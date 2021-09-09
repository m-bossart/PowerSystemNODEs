


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
        for (i,ω) in enumerate(V_freqs)
            val += V_coeffs[i][1]* sin.(ω * t)
            val += V_coeffs[i][2]* cos.(ω * t)
        end
        return val
    end
    θ_bias = get_internal_angle_bias(source)
    θ_freqs = get_internal_angle_frequencies(source)
    θ_coeffs = get_internal_angle_coefficients(source)
   function θ(t)
       val = θ_bias
       for (i,ω) in enumerate(θ_freqs)
           val += θ_coeffs[i][1]* sin.(ω * t)
           val += θ_coeffs[i][2]* cos.(ω * t)
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


function build_disturbances(sys)  #TODO make this more flexible, add options for which faults to include
    disturbances = []
    #BRANCH FAULTS
    lines = deepcopy(collect(get_components(Line, sys, x-> get_name(x) == "BUS 4       -BUS 5       -i_7"))) #TODO - change back to include all lines
    for l in lines
        push!(disturbances, BranchTrip(tfault,get_name(l)))
    end
    #REFERENCE CHANGE FAULTS
    injs = collect(get_components(DynamicInjection, sys,  x -> !(get_name(x) in get_name.(surrogate_gens))))
    for fault_inj in injs
        for Pref in Prefchange
            disturbance_ControlReferenceChange = ControlReferenceChange(tfault, fault_inj , PowerSimulationsDynamics.P_ref_index,  get_ext(fault_inj)["control_refs"][PowerSimulationsDynamics.P_ref_index] * Pref)
            #push!(disturbances, disturbance_ControlReferenceChange)
        end
    end
    return disturbances
end

function add_devices_to_surrogatize!(sys::System, n_devices::Integer, surrogate_bus_number::Integer, inf_bus_number:: Integer)
    param_range = (0.5, 2.0)
    surrogate_bus = collect(get_components(Bus,sys,x->get_number(x)==surrogate_bus_number))[1]
    inf_bus = collect(get_components(Bus,sys,x->get_number(x)==inf_bus_number))[1]

    surrogate_area = Area(;name = "surrogate")
    add_component!(sys,surrogate_area)
    set_area!(surrogate_bus, surrogate_area)
    set_area!(inf_bus, surrogate_area)

    gens = collect(get_components(ThermalStandard, sys, x->get_number(get_bus(x)) == surrogate_bus_number))

    !(length(gens) == 1) && @error "number of devices at surrogate bus not equal to one"
    gen = gens[1]
    total_rating = get_rating(gen) #doesn't impact dynamics
    total_base_power = get_base_power(gen)
    total_active_power = get_active_power(gen)
    remove_component!(sys,gen)
    for i in 1:n_devices
        g = ThermalStandard(
           name = string("gen",string(i)),
           available = true,
           status = true,
           bus = surrogate_bus,
           active_power = total_active_power, #Only divide base power by n_devices
           reactive_power = 0.0,
           rating =  total_rating/n_devices,
           active_power_limits=(min=0.0, max=3.0),
           reactive_power_limits= (-3.0,3.0),
           ramp_limits=nothing,
           operation_cost=ThreePartCost(nothing),
           base_power =  total_base_power/n_devices,
           )
       add_component!(sys, g)
       inv_typ = inv_case78(get_name(g))
       randomize_parameters!(inv_typ, param_range)
       add_component!(sys, inv_typ, g)
   end
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
       name = string("gen",string(1)),
       available = true,
       status = true,
       bus = collect(get_components(Bus,sys_init, x->get_bustype(x) == BusTypes.PV))[1],
       active_power = power_total / base_power_total, #Only divide base power by n_devices
       reactive_power = 0.0,
       rating =  base_power_total,
       active_power_limits = (0.0, 3.0),
       reactive_power_limits= (-3.0,3.0),
       ramp_limits=nothing,
       operation_cost=ThreePartCost(nothing),
       base_power =  base_power_total,
       )
    add_component!(sys_init, g)
    inv_typ = inv_case78(get_name(g))

    add_component!(sys_init, inv_typ, g)
    return sys_init
end

#NOTE The warning that the initialization fails in the source is because we just use the source to set the bus voltage.
#Doesn't make physical sense, but as long as the full system solves, it should be fine.
function initialize_sys!(sys::System, name::String, p)
    device = get_component(DynamicInverter, sys, name)
    set_parameters!(device, p)
    sim = Simulation!(
        MassMatrixModel,
        sys,
        pwd(),
        (0.0, 1.0),
    )
    x₀_dict = get_initial_conditions(sim)[get_name(device)]
    x₀ = [value for (key,value) in x₀_dict]
    refs = get_ext(device)["control_refs"]
    return x₀, refs
end

function set_bus_from_source(available_source::Source)
    Vsource = get_internal_voltage(available_source)
    set_magnitude!(get_bus(available_source),Vsource)
    θsource = get_internal_angle(available_source)
    set_angle!(get_bus(available_source),θsource)
end

function get_total_current_series(sim::Simulation)
    ir_total = []
    ii_total = []
    for (i,g) in enumerate(get_components(DynamicInjection, sys_train, x->typeof(x)!== PeriodicVariableSource))
        if i == 1
            ir_total = get_real_current_series(sim, get_name(g))[2]
            ii_total = get_imaginary_current_series(sim, get_name(g))[2]
        else
            ir_total .+= get_real_current_series(sim, get_name(g))[2]
            ii_total .+= get_imaginary_current_series(sim, get_name(g))[2]
        end
    end
    data_array =  zeros(Float64, (2, length(ir_total)))
    data_array[1,:] .= ir_total
    data_array[2,:] .= ii_total
    return data_array
end

function plot_pvs(tsteps, pvs::PeriodicVariableSource)
    V = zeros(length(tsteps))
    V = V .+ get_internal_voltage_bias(pvs)
    @info V
    retrieved_freqs = get_internal_voltage_frequencies(pvs)
    coeffs = get_internal_voltage_coefficients(pvs)
    for (i,ω) in enumerate(retrieved_freqs)
        V += coeffs[i][1]* sin.(ω.* tsteps)
        V += coeffs[i][2]* cos.(ω.* tsteps)
    end

    θ = zeros(length(tsteps))
    θ = θ .+ get_internal_angle_bias(pvs)
    @info θ
    retrieved_freqs = get_internal_angle_frequencies(pvs)
    coeffs = get_internal_angle_coefficients(pvs)
    for (i,ω) in enumerate(retrieved_freqs)
        θ += coeffs[i][1]* sin.(ω .* tsteps)
        θ += coeffs[i][2]* cos.(ω.* tsteps)
    end
    p1 = plot(tsteps,V, label = "plot from pvs coefficients")
    p2 = plot(tsteps,θ, label = "plot from pvs coefficients")
    return p1, p2
end


#Before you can initialize your surrogate, you need the true response in steady state.

#function Initialize(V,θ,Ir,Ii)

function cb_gfm_plot(sol)
    p1 = scatter(tsteps, sol[5,:], markersize = 2,  xaxis=:log,  label = "real current prediction")
    plot!(p1, ir_truth, label = "real current true")
    p2 = scatter(tsteps, sol[19,:], markersize = 2, xaxis=:log, label = "imag current prediction")
    plot!(p2, ii_truth, label = "imag current true")
    plt = plot(p1,p2,  layout=(2,1))
    push!(list_plots, plt)
    display_plots && display(plt)
end

function cb_gfm_nn_plot(pred, batch, time_batch)
    p1 = scatter(time_batch, pred[1,:], markersize=2, label = "real current prediction")
    plot!(p1, time_batch, batch[1,:],  xaxis=:log,  label = "real current true")
    p2 = scatter(time_batch, pred[2,:], markersize=2, label = "imag current prediction")
    plot!(p2, time_batch, batch[2,:],  xaxis=:log, label = "imag current true")
    plt = plot(p1,p2,layout=(2,1))
    push!(list_plots, plt)
    display_plots && display(plt)
end


function extending_ranges(datasize::Integer, groupsize::Integer)
    1 <= groupsize <= datasize || throw(
        DomainError(
            groupsize,
            "datasize must be positive and groupsize must to be within [1, datasize]",
        ),
    )
    return [1:min(datasize, i + groupsize - 1) for i in 1:groupsize:datasize]
end


function build_nn(input_dim, output_dim, nn_width, nn_hidden, nn_activation)
    if nn_hidden == 1 
        nn = FastChain(FastDense(input_dim, nn_width, nn_activation),
        FastDense(nn_width, nn_width, nn_activation),
        FastDense(nn_width, output_dim))
        return nn 
    elseif nn_hidden == 2 
        nn = FastChain(FastDense(input_dim, nn_width, nn_activation),
        FastDense(nn_width, nn_width, nn_activation),
        FastDense(nn_width, nn_width, nn_activation),
        FastDense(nn_width, output_dim))
        return nn 
    elseif nn_hidden == 3 
        nn = FastChain(FastDense(input_dim, nn_width, nn_activation),
        FastDense(nn_width, nn_width, nn_activation),
        FastDense(nn_width, nn_width, nn_activation),
        FastDense(nn_width, nn_width, nn_activation),
        FastDense(nn_width, output_dim))
        return nn 
    elseif nn_hidden == 4
        nn = FastChain(FastDense(input_dim, nn_width, nn_activation),
        FastDense(nn_width, nn_width, nn_activation),
        FastDense(nn_width, nn_width, nn_activation),
        FastDense(nn_width, nn_width, nn_activation),
        FastDense(nn_width, nn_width, nn_activation),
        FastDense(nn_width, output_dim))
        return nn 
    elseif nn_hidden == 5
        nn = FastChain(FastDense(input_dim, nn_width, nn_activation),
        FastDense(nn_width, nn_width, nn_activation),
        FastDense(nn_width, nn_width, nn_activation),
        FastDense(nn_width, nn_width, nn_activation),
        FastDense(nn_width, nn_width, nn_activation),
        FastDense(nn_width, nn_width, nn_activation),
        FastDense(nn_width, output_dim))
        return nn 
    else 
        @error "build_nn does not support the provided nn depth"
        return false 
    end 
end 