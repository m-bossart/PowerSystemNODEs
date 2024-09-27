using PowerSystems
using Plots

##
# Option 1: Uniform double of line impedances in the surrogate portion of the system. 
sys = System(joinpath(pwd(), "systems", "36bus_fix.json"))
surrogate_buses = vcat(21:29, 31:39)
for b in get_components(DynamicBranch, sys)
    i = b.branch
    to_bus_number = get_number(get_to(get_arc(i)))
    from_bus_number =  get_number(get_from(get_arc(i)))
    if (to_bus_number in surrogate_buses) && (from_bus_number in surrogate_buses) 
        @show from_bus_number, to_bus_number

        set_x!(i, get_x(i) * 4)
        set_r!(i, get_r(i) * 4)
    end 
end 
to_json(sys, joinpath(pwd(), "systems", "36bus_4x_line_Z.json"))



##
# Option 2: Change parameters: Chance 
sys = System(joinpath(pwd(), "systems", "36bus_fix.json"))

for i in get_components(DynamicInverter{AverageConverter, OuterControl{ActivePowerDroop, ReactivePowerDroop}, VoltageModeControl, FixedDCSource, FixedFrequency, LCLFilter}, sys)    
    oc = get_outer_control(i)
    a = PowerSystems.get_active_power_control(oc)
    set_Rp!(a, get_Rp(a) * rand(Distributions.Uniform(0.33, 3.0)))
    set_ωz!(a, get_ωz(a) *  rand(Distributions.Uniform(0.33, 3.0)))

    r = PowerSystems.get_reactive_power_control(oc)
    set_kq!(r, get_kq(r) * rand(Distributions.Uniform(0.33, 3.0)))
    set_ωf!(r, get_ωf(r) * rand(Distributions.Uniform(0.33, 3.0)))
end 

for i in get_components(DynamicInverter{AverageConverter, OuterControl{ActivePowerPI, ReactivePowerPI}, CurrentModeControl, FixedDCSource, KauraPLL, LCLFilter}, sys)
    oc = get_outer_control(i)
    a = PowerSystems.get_active_power_control(oc)
    set_Kp_p!(a, get_Kp_p(a) *  rand(Distributions.Uniform(0.33, 3.0)))
    set_Ki_p!(a, get_Ki_p(a) * rand(Distributions.Uniform(0.33, 3.0)))
    set_ωz!(a, get_ωz(a) *  rand(Distributions.Uniform(0.33, 3.0)))

    r = PowerSystems.get_reactive_power_control(oc)
    set_Kp_q!(r, get_Kp_q(r) * rand(Distributions.Uniform(0.33, 3.0)))
    set_Ki_q!(r, get_Ki_q(r) *  rand(Distributions.Uniform(0.33, 3.0)))
    set_ωf!(r, get_ωf(r) *  rand(Distributions.Uniform(0.33, 3.0)))
end 
to_json(sys, joinpath(pwd(), "systems", "36bus_increase_parameter_variation.json"))



