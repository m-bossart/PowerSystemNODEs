using PowerSystems
const PSY = PowerSystems

#CLASSICAL GENERATOR CONSTRUCTION
machine_classic() = BaseMachine(0.0, 0.2995, 0.7087) #ra, Xd_p, eq_p (eq_p will change)
shaft_damping() = SingleMass(3.148, 2.0)
avr_none() = AVRFixed(0.0)
tg_none() = TGFixed(1.0) #efficiency
pss_none() = PSSFixed(0.0)
function dyn_gen_classic(generator)    #1.0 is ωref
    return DynamicGenerator(
        generator,
        1.0,
        machine_classic(),
        shaft_damping(),
        avr_none(),
        tg_none(),
        pss_none(),
    )
end

#################################### CONVERTER #########################################################
converter_low_power() = AverageConverter(rated_voltage = 690.0, rated_current = 2.75)
converter_high_power() = AverageConverter(rated_voltage = 138.0, rated_current = 100.0)

#################################### OUTER CONTROL #####################################################
function outer_control()
    function virtual_inertia()
        return VirtualInertia(Ta = 2.0, kd = 400.0, kω = 20.0)
    end
    function reactive_droop()
        return ReactivePowerDroop(kq = 0.2, ωf = 1000.0)
    end
    return OuterControl(virtual_inertia(), reactive_droop())
end

function outer_control_droop()
    function active_droop()
        return PSY.ActivePowerDroop(Rp = 0.05, ωz = 2 * pi * 5)
    end
    function reactive_droop()
        return ReactivePowerDroop(kq = 0.2, ωf = 1000.0)
    end
    return OuterControl(active_droop(), reactive_droop())
end

function outer_control_gfoll()
    function active_pi()
        return ActivePowerPI(Kp_p = 2.0, Ki_p = 30.0, ωz = 0.132 * 2 * pi * 50)
    end
    function reactive_pi()
        return ReactivePowerPI(Kp_q = 2.0, Ki_q = 30.0, ωf = 0.132 * 2 * pi * 50.0)
    end
    return OuterControl(active_pi(), reactive_pi())
end

#################################### INNER CONRTOL ######################################################
inner_control() = VoltageModeControl(
    kpv = 0.59,     #Voltage controller proportional gain
    kiv = 736.0,    #Voltage controller integral gain
    kffv = 0.0,     #Binary variable enabling the voltage feed-forward in output of current controllers
    rv = 0.0,       #Virtual resistance in pu
    lv = 0.2,       #Virtual inductance in pu
    kpc = 1.27,     #Current controller proportional gain
    kic = 14.3,     #Current controller integral gain
    kffi = 0.0,     #Binary variable enabling the current feed-forward in output of current controllers
    ωad = 50.0,     #Active damping low pass filter cut-off frequency
    kad = 0.2,      #Active damping gain
)
current_mode_inner() = CurrentModeControl(
    kpc = 0.37,     #Current controller proportional gain
    kic = 0.7,     #Current controller integral gain
    kffv = 1.0,     #Binary variable enabling the voltage feed-forward in output of current controllers
)

#################################### DC SOURCE #########################################################
dc_source_lv() = FixedDCSource(voltage = 600.0)

#################################### FREQUENCY ESTIMATOR ###############################################
pll() = KauraPLL(
    ω_lp = 500.0, #Cut-off frequency for LowPass filter of PLL filter.
    kp_pll = 0.084,  #PLL proportional gain
    ki_pll = 4.69,   #PLL integral gain
)
reduced_pll() = ReducedOrderPLL(
    ω_lp = 1.32 * 2 * pi * 50, #Cut-off frequency for LowPass filter of PLL filter.
    kp_pll = 2.0,  #PLL proportional gain
    ki_pll = 20.0,   #PLL integral gain
)
no_pll() = PSY.FixedFrequency()

#################################### LCL FILTERS ######################################################
filt() = LCLFilter(lf = 0.08, rf = 0.003, cf = 0.074, lg = 0.2, rg = 0.01)
filt_gfoll() = LCLFilter(lf = 0.009, rf = 0.016, cf = 2.5, lg = 0.002, rg = 0.003)

#################################### FULL INVERTERS ###################################################

function inv_case78(static_device)
    return DynamicInverter(
        static_device,
        1.0, # ω_ref,
        converter_high_power(), #converter
        outer_control(), #outer control
        inner_control(), #inner control voltage source
        dc_source_lv(), #dc source
        pll(), #pll
        filt(), #filter
    )
end
function inv_darco_droop(static_device)
    return DynamicInverter(
        static_device,
        1.0, #ω_ref
        converter_low_power(), #converter
        outer_control_droop(), #outercontrol
        inner_control(), #inner_control
        dc_source_lv(),
        no_pll(),
        filt(),
    ) #pss
end
function inv_gfoll(static_device)
    return DynamicInverter(
        static_device,
        1.0, #ω_ref
        converter_low_power(), #converter
        outer_control_gfoll(), #outercontrol
        current_mode_inner(), #inner_control
        dc_source_lv(),
        reduced_pll(),
        filt_gfoll(),
    ) #pss
end

function add_source_to_ref(sys::System)
    for g in get_components(StaticInjection, sys)
        isa(g, ElectricLoad) && continue
        g.bus.bustype == BusTypes.REF &&
            error("A device is already attached to the REF bus")
    end
    slack_bus = [b for b in get_components(Bus, sys) if b.bustype == BusTypes.REF][1]
    inf_source = Source(
        name = "InfBus", #name
        active_power = 0.0,
        available = true, #availability
        reactive_power = 0.0,
        bus = slack_bus, #bus
        R_th = 0.0, #Rth
        X_th = 5e-6, #Xth
    )
    add_component!(sys, inf_source)
    return
end
