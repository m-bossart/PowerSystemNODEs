
"""
    randomize_inv_parameters!(inv::DynamicInverter, param_range::Tuple{Float64,Float64})

Change each inverter parameter by scaling by randomly generated factor in param_range
"""
function randomize_inv_parameters!(inv::DynamicInverter, param_range::Tuple{Float64,Float64})
    #pll
    set_ω_lp!(inv.freq_estimator, get_ω_lp(inv.freq_estimator)*rand(Uniform(param_range[1],param_range[2])))
    set_kp_pll!(inv.freq_estimator, get_kp_pll(inv.freq_estimator)*rand(Uniform(param_range[1],param_range[2])))
    set_ki_pll!(inv.freq_estimator, get_ki_pll(inv.freq_estimator)*rand(Uniform(param_range[1],param_range[2])))
    #outer control
    PSY.set_Ta!(inv.outer_control.active_power, PSY.get_Ta(inv.outer_control.active_power)*rand(Uniform(param_range[1],param_range[2])))
    PSY.set_kd!(inv.outer_control.active_power, PSY.get_kd(inv.outer_control.active_power)*rand(Uniform(param_range[1],param_range[2])))
    PSY.set_kω!(inv.outer_control.active_power, PSY.get_kω(inv.outer_control.active_power)*rand(Uniform(param_range[1],param_range[2])))
    PSY.set_kq!(inv.outer_control.reactive_power, PSY.get_kq(inv.outer_control.reactive_power)*rand(Uniform(param_range[1],param_range[2])))
    PSY.set_ωf!(inv.outer_control.reactive_power, PSY.get_ωf(inv.outer_control.reactive_power)*rand(Uniform(param_range[1],param_range[2])))
    #inner control
    set_kpv!(inv.inner_control, get_kpv(inv.inner_control)*rand(Uniform(param_range[1],param_range[2])))
    set_kiv!(inv.inner_control, get_kiv(inv.inner_control)*rand(Uniform(param_range[1],param_range[2])))
    set_kffv!(inv.inner_control, get_kffv(inv.inner_control)*rand(Uniform(param_range[1],param_range[2])))
    set_rv!(inv.inner_control, get_rv(inv.inner_control)*rand(Uniform(param_range[1],param_range[2])))
    set_lv!(inv.inner_control, get_lv(inv.inner_control)*rand(Uniform(param_range[1],param_range[2])))
    set_kpc!(inv.inner_control, get_kpc(inv.inner_control)*rand(Uniform(param_range[1],param_range[2])))
    set_kic!(inv.inner_control, get_kic(inv.inner_control)*rand(Uniform(param_range[1],param_range[2])))
    set_kffi!(inv.inner_control, get_kffi(inv.inner_control)*rand(Uniform(param_range[1],param_range[2])))
    set_ωad!(inv.inner_control, get_ωad(inv.inner_control)*rand(Uniform(param_range[1],param_range[2])))
    set_kad!(inv.inner_control, get_kad(inv.inner_control)*rand(Uniform(param_range[1],param_range[2])))
    #lcl
    set_lf!(inv.filter, get_lf(inv.filter)*rand(Uniform(param_range[1],param_range[2])))
    set_rf!(inv.filter, get_rf(inv.filter)*rand(Uniform(param_range[1],param_range[2])))
    set_cf!(inv.filter, get_cf(inv.filter)*rand(Uniform(param_range[1],param_range[2])))
    set_lg!(inv.filter, get_lg(inv.filter)*rand(Uniform(param_range[1],param_range[2])))
    set_rg!(inv.filter, get_rg(inv.filter)*rand(Uniform(param_range[1],param_range[2])))
end


function set_inv_parameters!(inv::DynamicInverter, p::Vector{Float64})
    #pll
    set_ω_lp!(inv.freq_estimator, p[1])
    set_kp_pll!(inv.freq_estimator, p[2])
    set_ki_pll!(inv.freq_estimator ,p[3])
    #outer control
    PSY.set_Ta!(inv.outer_control.active_power, p[4])
    PSY.set_kd!(inv.outer_control.active_power, p[5])
    PSY.set_kω!(inv.outer_control.active_power, p[6])
    PSY.set_kq!(inv.outer_control.reactive_power, p[7])
    PSY.set_ωf!(inv.outer_control.reactive_power,p[8])
    #inner control
    set_kpv!(inv.inner_control, p[9])
    set_kiv!(inv.inner_control, p[10])
    set_kffv!(inv.inner_control, p[11])
    set_rv!(inv.inner_control, p[12])
    set_lv!(inv.inner_control, p[13])
    set_kpc!(inv.inner_control, p[14])
    set_kic!(inv.inner_control, p[15])
    set_kffi!(inv.inner_control, p[16])
    set_ωad!(inv.inner_control, p[17])
    set_kad!(inv.inner_control, p[18])
    #lcl
    set_lf!(inv.filter, p[19])
    set_rf!(inv.filter, p[20])
    set_cf!(inv.filter, p[21])
    set_lg!(inv.filter, p[22])
    set_rg!(inv.filter, p[23])
end
"""
    build_train_test(sys_faults::System, sys_structure::System, train_split)

Builds a train and test system by combining a system with pre-defined faults and a system with the structure
"""

function build_train_test(sys_faults::System,  sys_full::System, Ref_bus_number::Integer, train_split; add_pvs = true )
    sys_train = deepcopy(sys_full)
    #remove_components!(sys_train, FixedAdmittance) #BUG add back if you include fixed admittance
    remove_components!(sys_train, PowerLoad)
    remove_components!(sys_train, LoadZone)
    remove_components!( x-> !(get_name(get_area(get_to(x))) == "surrogate" && get_name(get_area(get_from(x)))  == "surrogate"), sys_train, Arc)
    remove_components!( x-> !(get_name(get_area(get_to(get_arc(x)))) == "surrogate" && get_name(get_area(get_from(get_arc(x))))  == "surrogate"), sys_train, Transformer2W)
    remove_components!( x-> !(get_name(get_area(get_to(get_arc(x)))) == "surrogate" && get_name(get_area(get_from(get_arc(x))))  == "surrogate"), sys_train, Line)
    gens_to_remove = get_components(ThermalStandard,sys_train,  x-> !(get_name(get_area(get_bus(x))) == "surrogate"))
    for g in gens_to_remove
        dyn = get_dynamic_injector(g)
        (dyn !== nothing ) && remove_component!(sys_train, dyn)
        remove_component!(sys_train, g)
    end
    remove_components!( x-> !(get_name(get_area((x))) == "surrogate"), sys_train, Bus)
    @info length(collect(get_components(Bus, sys_train)))
    #remove_components!(sys_train, FixedAdmittance)

    #Remove all buses and
    slack_bus_train = collect(get_components(Bus, sys_train, x-> get_number(x) == Ref_bus_number))[1]
    set_bustype!(slack_bus_train,BusTypes.REF)

    sys_test = deepcopy(sys_train)
    slack_bus_test = collect(get_components(Bus,sys_test, x->get_bustype(x) == BusTypes.REF))[1]

    sources = get_components(Source, sys_faults)

    for s in sources
        pvs = get_dynamic_injector(s)
        remove_component!(sys_faults, pvs)
        remove_component!(sys_faults, s)

        if rand()<train_split
            set_bus!(s,slack_bus_train)
            display(s)
            add_component!(sys_train, s)
            add_pvs && add_component!(sys_train, pvs, s)
        else
            set_bus!(s,slack_bus_test)
            add_component!(sys_test, s)
            add_pvs && add_component!(sys_test, pvs, s)
        end
    end
    return sys_train, sys_test
end

"""
    PVS_to_function_of_time(source::PeriodicVariableSource)

Takes in a PeriodicVariableSource from PowerSystems and generates functions of time for voltage magnitude and angle
"""
function Source_to_function_of_time(source::PeriodicVariableSource)
     V_bias = get_internal_voltage_bias(source)
     V_freqs = get_internal_voltage_frequencies(source)
     V_coeffs = get_internal_voltage_coefficients(source)
    function V(t)
        val = V_bias
        for (i,f) in enumerate(V_freqs)
            val += V_coeffs[i][1]* sin.(f *2 * pi * t)
            val += V_coeffs[i][2]* cos.(f *2 * pi * t)
        end
        return val
    end
    θ_bias = get_internal_angle_bias(source)
    θ_freqs = get_internal_angle_frequencies(source)
    θ_coeffs = get_internal_angle_coefficients(source)
   function θ(t)
       val = θ_bias
       for (i,f) in enumerate(θ_freqs)
           val += θ_coeffs[i][1]* sin.(f *2 * pi * t)
           val += θ_coeffs[i][2]* cos.(f *2 * pi * t)
       end
       return val
   end
    return (V, θ)
end

function Source_to_function_of_time(source::Source)
    function V(t)
        return get_internal_voltage(source)
    end

   function θ(t)
       return get_internal_angle(source)
   end
    return (V, θ)
end




"""
    activate_next_source!(sys::System)

Either activate the first source if none are available, or make the next source available.
To be used in training surrogate to move on to the next system disturbance. Returns the available source
"""
function activate_next_source!(sys::System)
    all_sources = collect(get_components(Source,sys))
    active_sources = collect(get_components(Source,sys, x -> PSY.get_available(x)))
    if length(active_sources) < 1
        @info "no active sources in the system, activating the first source"
        first_source = collect(get_components(Source,sys))[1]
        set_available!(first_source, true)
        return first_source
    elseif length(active_sources) > 1
        @error "more than one active source, cannot determine next active source"
    else
        for (i,source) in enumerate(all_sources)
            if active_sources[1] == source
                set_available!(all_sources[i], false)
                if source !== last(all_sources)
                    set_available!(all_sources[i+1], true )
                    @info "found active source, setting next source active"
                    return all_sources[i+1]
                else
                    set_available!(all_sources[1], true)
                    @info "the last source is active, starting over at index 1 "
                    return all_sources[1]
                end
            end
        end
    end
end

"""

Makes a Float64 Mass Matrix of ones for the ODEProblem. Takes # of differential and algebraic states

"""
function MassMatrix(n_differential::Integer, n_algebraic::Integer)
    n_states = n_differential + n_algebraic
    M = Float64.(zeros(n_states,n_states))
    for i = 1:n_differential   #-2     Include if using the IB version (last two equations are algebraic)
      M[i,i] = 1.0
    end
    return M
end


function find_acbranch(from_bus_number::Int, to_bus_number::Int)
    for b in get_components(ACBranch,sys)
        (b.arc.from.number == from_bus_number) && ( b.arc.to.number == to_bus_number) && return(b)
    end
end


function build_sys_init(sys_train::System)
    sys_init = deepcopy(sys_train)
    base_power_total = 0.0
    power_total = 0.0
    for gfm in get_components(ThermalStandard,sys_init, x->typeof(get_dynamic_injector(x)) == DynamicInverter{AverageConverter, OuterControl{VirtualInertia, ReactivePowerDroop}, VoltageModeControl, FixedDCSource, KauraPLL, LCLFilter})
        base_power_total += get_base_power(gfm)
        power_total +=  get_base_power(gfm) * get_active_power(gfm)
        @info base_power_total
        @info power_total
        remove_component!(sys_init, get_dynamic_injector(gfm))
        remove_component!(sys_init, gfm)
    end
    g = ThermalStandard(
       name = string(1),
       available = true,
       status = true,
       bus = collect(get_components(Bus,sys_init, x->get_bustype(x) == BusTypes.PV))[1],
       active_power = power_total / base_power_total, #Only divide base power by n_devices
       reactive_power = 0.0,
       rating =  total_rating/n_devices,
       active_power_limits=(min=0.0, max=3.0),
       reactive_power_limits= (-3.0,3.0),
       ramp_limits=nothing,
       operation_cost=ThreePartCost(nothing),
       base_power =  base_power_total,
       )
    add_component!(sys_init, g)
    inv_typ = inv_case78(get_name(g))
    #set_inv_parameters!(inv_typ, p_start)
    add_component!(sys_init, inv_typ, g)
    return sys_init
end

function initalize_sys_init!(sys::System, p::Vector{Float64})
    gfm = collect(get_components(DynamicInverter, sys))[1]
    #@info "# of gfms", get_components(DynamicInverter, sys)
    set_inv_parameters!(gfm, p)
    sim = Simulation!(
        MassMatrixModel,
        sys,
        pwd(),
        (0.0, 1.0),
    )
    x₀_dict = get_initial_conditions(sim)[get_name(gfm)]
    x₀ = Float64.([value for (key,value) in x₀_dict])
    refs = get_ext(gfm)["control_refs"]
    return x₀, refs
end

function set_bus_from_source(available_source::Source)
    Vsource = get_internal_voltage(available_source)
    set_magnitude!(get_bus(available_source),Vsource)
    θsource = get_internal_angle(available_source)
    set_angle!(get_bus(available_source),θsource)
end


#Before you can initialize your surrogate, you need the true response in steady state.

#function Initialize(V,θ,Ir,Ii)
