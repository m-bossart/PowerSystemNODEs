
"""
    build_train_test(sys_faults::System, sys_structure::System)

Builds a train system by combining a system with pre-defined faults and a system with the structure
"""

function build_sys_train(sys_faults::System, sys_full::System, Ref_bus_number::Integer)
    sys_train = deepcopy(sys_full)
    #remove_components!(sys_train, FixedAdmittance) #BUG add back if you include fixed admittance
    remove_components!(sys_train, PowerLoad)
    remove_components!(sys_train, LoadZone)
    remove_components!(
        x ->
            !(
                get_name(get_area(get_to(x))) == "surrogate" &&
                get_name(get_area(get_from(x))) == "surrogate"
            ),
        sys_train,
        Arc,
    )
    remove_components!(
        x ->
            !(
                get_name(get_area(get_to(get_arc(x)))) == "surrogate" &&
                get_name(get_area(get_from(get_arc(x)))) == "surrogate"
            ),
        sys_train,
        Transformer2W,
    )
    remove_components!(
        x ->
            !(
                get_name(get_area(get_to(get_arc(x)))) == "surrogate" &&
                get_name(get_area(get_from(get_arc(x)))) == "surrogate"
            ),
        sys_train,
        Line,
    )
    gens_to_remove = get_components(
        ThermalStandard,
        sys_train,
        x -> !(get_name(get_area(get_bus(x))) == "surrogate"),
    )
    for g in gens_to_remove
        dyn = get_dynamic_injector(g)
        (dyn !== nothing) && remove_component!(sys_train, dyn)
        remove_component!(sys_train, g)
    end
    remove_components!(x -> !(get_name(get_area((x))) == "surrogate"), sys_train, Bus)
    @info length(collect(get_components(Bus, sys_train)))
    #remove_components!(sys_train, FixedAdmittance)

    #Remove all buses and
    slack_bus_train =
        collect(get_components(Bus, sys_train, x -> get_number(x) == Ref_bus_number))[1]
    set_bustype!(slack_bus_train, BusTypes.REF)

    #sys_test = deepcopy(sys_train)
    #slack_bus_test =
    #   collect(get_components(Bus, sys_test, x -> get_bustype(x) == BusTypes.REF))[1]

    sources = get_components(Source, sys_faults)

    for s in sources
        pvs = get_dynamic_injector(s)
        remove_component!(sys_faults, pvs)
        remove_component!(sys_faults, s)

        set_bus!(s, slack_bus_train)
        display(s)
        add_component!(sys_train, s)
        add_component!(sys_train, pvs, s)
    end
    return sys_train
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
        for (i, ω) in enumerate(V_freqs)
            val += V_coeffs[i][1] * sin.(ω * t)
            val += V_coeffs[i][2] * cos.(ω * t)
        end
        return val
    end
    θ_bias = get_internal_angle_bias(source)
    θ_freqs = get_internal_angle_frequencies(source)
    θ_coeffs = get_internal_angle_coefficients(source)
    function θ(t)
        val = θ_bias
        for (i, ω) in enumerate(θ_freqs)
            val += θ_coeffs[i][1] * sin.(ω * t)
            val += θ_coeffs[i][2] * cos.(ω * t)
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

function pad_tanh(t, y)
    first = y[1]
    last = y[end]
    Δt = t[2] - t[1] #assume linear spacing
    n = length(t)
    A = (y[1] - y[end]) / 2
    C = y[end] + A

    B = 20 / (t[end] - t[1])
    @info A, B, C
    @info y[end]
    @info y[1]
    t_add = (t[end] + Δt):Δt:(t[end] + (t[end] - t[1]))
    y_add = A * tanh.((t_add .- t_add[Int(length(t_add) / 2)]) .* B) .+ C
    t_new = vcat(t, t_add)
    y_new = vcat(y, y_add)
    return (t_new, y_new)
end

function build_disturbances(sys)  #TODO make this more flexible, add options for which faults to include
    disturbances = []
    #BRANCH FAULTS
    lines = deepcopy(
        collect(get_components(Line, sys, x -> get_name(x) == "BUS 13-BUS 14-i_17")),
    ) #"BUS 4-BUS 5-i_7"
    println(lines)
    for l in lines
        #push!(disturbances, BranchTrip(tfault,Line,get_name(l)))
    end
    #REFERENCE CHANGE FAULTS
    injs = collect(
        get_components(
            DynamicInjection,
            sys,
            x -> !(get_name(x) in get_name.(surrogate_gens)),
        ),
    )
    for fault_inj in injs
        for Pref in Prefchange
            disturbance_ControlReferenceChange = ControlReferenceChange(
                tfault,
                fault_inj,
                :P_ref,
                get_P_ref(fault_inj) * Pref,
            )

            (get_name(fault_inj) == "generator-15-1") &&
                push!(disturbances, disturbance_ControlReferenceChange)
        end
    end
    return disturbances
end

function add_devices_to_surrogatize!(
    sys::System,
    n_devices::Integer,
    surrogate_bus_number::Integer,
    inf_bus_number::Integer,
)
    param_range = (0.5, 2.0)
    surrogate_bus =
        collect(get_components(Bus, sys, x -> get_number(x) == surrogate_bus_number))[1]
    inf_bus = collect(get_components(Bus, sys, x -> get_number(x) == inf_bus_number))[1]

    surrogate_area = Area(; name = "surrogate")
    add_component!(sys, surrogate_area)
    set_area!(surrogate_bus, surrogate_area)
    set_area!(inf_bus, surrogate_area)

    gens = collect(
        get_components(
            ThermalStandard,
            sys,
            x -> get_number(get_bus(x)) == surrogate_bus_number,
        ),
    )

    !(length(gens) == 1) && @error "number of devices at surrogate bus not equal to one"
    gen = gens[1]
    total_rating = get_rating(gen) #doesn't impact dynamics
    total_base_power = get_base_power(gen)
    total_active_power = get_active_power(gen)
    remove_component!(sys, gen)
    for i in 1:n_devices
        g = ThermalStandard(
            name = string("gen", string(i)),
            available = true,
            status = true,
            bus = surrogate_bus,
            active_power = total_active_power, #Only divide base power by n_devices
            reactive_power = 0.0,
            rating = total_rating / n_devices,
            active_power_limits = (min = 0.0, max = 3.0),
            reactive_power_limits = (-3.0, 3.0),
            ramp_limits = nothing,
            operation_cost = ThreePartCost(nothing),
            base_power = total_base_power / n_devices,
        )
        add_component!(sys, g)
        if (i == 1)
            inv_typ = inv_case78(get_name(g))
            add_component!(sys, inv_typ, g)
        end
        if (i == 2)
            inv_typ = inv_gfoll(get_name(g))
            add_component!(sys, inv_typ, g)
        end
        if (i == 3)
            inv_typ = inv_darco_droop(get_name(g))
            add_component!(sys, inv_typ, g)
        end
    end
end

"""
    activate_next_source!(sys::System)

Either activate the first source if none are available, or make the next source available.
To be used in training surrogate to move on to the next system disturbance. Returns the available source
"""
function activate_next_source!(sys::System)
    all_sources = collect(get_components(Source, sys))
    active_sources = collect(get_components(Source, sys, x -> PSY.get_available(x)))
    if length(active_sources) < 1
        @info "no active sources in the system, activating the first source"
        first_source = collect(get_components(Source, sys))[1]
        set_available!(first_source, true)
        return first_source
    elseif length(active_sources) > 1
        @error "more than one active source, cannot determine next active source"
    else
        for (i, source) in enumerate(all_sources)
            if active_sources[1] == source
                set_available!(all_sources[i], false)
                if source !== last(all_sources)
                    set_available!(all_sources[i + 1], true)
                    @info "found active source, setting next source active"
                    return all_sources[i + 1]
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
    M = Float64.(zeros(n_states, n_states))
    for i in 1:n_differential
        M[i, i] = 1.0
    end
    return M
end

function find_acbranch(from_bus_number::Int, to_bus_number::Int)
    for b in get_components(ACBranch, sys)
        (b.arc.from.number == from_bus_number) &&
            (b.arc.to.number == to_bus_number) &&
            return (b)
    end
end

"""
Test function description
"""
function build_sys_init(sys_train::System)
    sys_init = deepcopy(sys_train)
    base_power_total = 0.0
    power_total = 0.0
    gfms = collect(get_components(ThermalStandard, sys_init))

    for gfm in gfms
        base_power_total += get_base_power(gfm)
        power_total += get_base_power(gfm) * get_active_power(gfm)
        @info base_power_total
        @info power_total
        remove_component!(sys_init, get_dynamic_injector(gfm))
        remove_component!(sys_init, gfm)
    end
    g = ThermalStandard(
        name = string("gen", string(1)),
        available = true,
        status = true,
        bus = collect(get_components(Bus, sys_init, x -> get_bustype(x) == BusTypes.PV))[1],
        active_power = power_total / base_power_total, #Only divide base power by n_devices
        reactive_power = 0.0,
        rating = base_power_total,
        active_power_limits = (0.0, 3.0),
        reactive_power_limits = (-3.0, 3.0),
        ramp_limits = nothing,
        operation_cost = ThreePartCost(nothing),
        base_power = base_power_total,
    )
    add_component!(sys_init, g)
    inv_typ = inv_case78(get_name(g))
    add_component!(sys_init, inv_typ, g)
    p_inv = get_parameters(inv_typ)
    return sys_init, p_inv
end

function get_wrapper(wrappers, name)
    for w in wrappers
        if PSY.get_name(w) == name
            return w
        end
    end
end
#NOTE The warning that the initialization fails in the source is because we just use the source to set the bus voltage.
#Doesn't make physical sense, but as long as the full system solves, it should be fine.
function initialize_sys!(sys::System, name::String)
    device = get_component(DynamicInverter, sys, name)
    bus = get_bus(get_component(StaticInjection, sys, name)).number
    # set_parameters!(device, p)
    sim = Simulation!(
        MassMatrixModel,
        sys,
        pwd(),
        (0.0, 1.0),
        console_level = PSID_CONSOLE_LEVEL,
        file_level = PSID_FILE_LEVEL,
    )
    x₀_dict = read_initial_conditions(sim)[get_name(device)]
    x₀ = [value for (key, value) in x₀_dict]
    wrappers = sim.inputs.dynamic_injectors
    w = get_wrapper(wrappers, "gen1")
    refs = [w.V_ref.x, w.ω_ref.x, w.P_ref.x, w.Q_ref.x]
    Vr0 = read_initial_conditions(sim)["V_R"][bus]
    Vi0 = read_initial_conditions(sim)["V_I"][bus]
    return x₀, refs, Vr0, Vi0
end

function set_bus_from_source(available_source::Source)
    Vsource = get_internal_voltage(available_source)
    set_magnitude!(get_bus(available_source), Vsource)
    θsource = get_internal_angle(available_source)
    set_angle!(get_bus(available_source), θsource)
end

function get_total_current_series(sim)
    ir_total = []
    ii_total = []
    for (i, g) in enumerate(
        get_components(
            DynamicInjection,
            sim.sys,
            x -> typeof(x) !== PeriodicVariableSource,
        ),
    )
        results = read_results(sim)
        if i == 1
            ir_total = get_real_current_series(results, get_name(g))[2]
            ii_total = get_imaginary_current_series(results, get_name(g))[2]
        else
            ir_total .+= get_real_current_series(results, get_name(g))[2]
            ii_total .+= get_imaginary_current_series(results, get_name(g))[2]
        end
    end
    data_array = zeros(Float64, (2, length(ir_total)))
    data_array[1, :] .= ir_total
    data_array[2, :] .= ii_total
    return data_array
end

function plot_pvs(tsteps, pvs::PeriodicVariableSource, xaxis)
    V = zeros(length(tsteps))
    V = V .+ get_internal_voltage_bias(pvs)
    retrieved_freqs = get_internal_voltage_frequencies(pvs)
    coeffs = get_internal_voltage_coefficients(pvs)
    for (i, ω) in enumerate(retrieved_freqs)
        V += coeffs[i][1] * sin.(ω .* tsteps)
        V += coeffs[i][2] * cos.(ω .* tsteps)
    end

    θ = zeros(length(tsteps))
    θ = θ .+ get_internal_angle_bias(pvs)
    retrieved_freqs = get_internal_angle_frequencies(pvs)
    coeffs = get_internal_angle_coefficients(pvs)
    for (i, ω) in enumerate(retrieved_freqs)
        θ += coeffs[i][1] * sin.(ω .* tsteps)
        θ += coeffs[i][2] * cos.(ω .* tsteps)
    end
    p1 =
        plot(tsteps, V, label = "plot from pvs coefficients", xaxis = xaxis, legend = false)
    p2 = plot(tsteps, θ, label = "plot from pvs coefficients", xaxis = xaxis)
    return p1, p2
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
        nn = FastChain(
            FastDense(input_dim, nn_width, nn_activation),
            FastDense(nn_width, nn_width, nn_activation),
            FastDense(nn_width, output_dim),
        )
        return nn
    elseif nn_hidden == 2
        nn = FastChain(
            FastDense(input_dim, nn_width, nn_activation),
            FastDense(nn_width, nn_width, nn_activation),
            FastDense(nn_width, nn_width, nn_activation),
            FastDense(nn_width, output_dim),
        )
        return nn
    elseif nn_hidden == 3
        nn = FastChain(
            FastDense(input_dim, nn_width, nn_activation),
            FastDense(nn_width, nn_width, nn_activation),
            FastDense(nn_width, nn_width, nn_activation),
            FastDense(nn_width, nn_width, nn_activation),
            FastDense(nn_width, output_dim),
        )
        return nn
    elseif nn_hidden == 4
        nn = FastChain(
            FastDense(input_dim, nn_width, nn_activation),
            FastDense(nn_width, nn_width, nn_activation),
            FastDense(nn_width, nn_width, nn_activation),
            FastDense(nn_width, nn_width, nn_activation),
            FastDense(nn_width, nn_width, nn_activation),
            FastDense(nn_width, output_dim),
        )
        return nn
    elseif nn_hidden == 5
        nn = FastChain(
            FastDense(input_dim, nn_width, nn_activation),
            FastDense(nn_width, nn_width, nn_activation),
            FastDense(nn_width, nn_width, nn_activation),
            FastDense(nn_width, nn_width, nn_activation),
            FastDense(nn_width, nn_width, nn_activation),
            FastDense(nn_width, nn_width, nn_activation),
            FastDense(nn_width, output_dim),
        )
        return nn
    else
        @error "build_nn does not support the provided nn depth"
        return false
    end
end

function get_parameters(inv::DynamicInverter)
    #pll
    p = Vector{Float64}(undef, 23)
    p[1] = get_ω_lp(inv.freq_estimator)
    p[2] = get_kp_pll(inv.freq_estimator)
    p[3] = get_ki_pll(inv.freq_estimator)
    #outer control
    p[4] = PSY.get_Ta(inv.outer_control.active_power)
    p[5] = PSY.get_kd(inv.outer_control.active_power)
    p[6] = PSY.get_kω(inv.outer_control.active_power)
    p[7] = PSY.get_kq(inv.outer_control.reactive_power)
    p[8] = PSY.get_ωf(inv.outer_control.reactive_power)
    #inner control
    p[9] = get_kpv(inv.inner_control)
    p[10] = get_kiv(inv.inner_control)
    p[11] = get_kffv(inv.inner_control)
    p[12] = get_rv(inv.inner_control)
    p[13] = get_lv(inv.inner_control)
    p[14] = get_kpc(inv.inner_control)
    p[15] = get_kic(inv.inner_control)
    p[16] = get_kffi(inv.inner_control)
    p[17] = get_ωad(inv.inner_control)
    p[18] = get_kad(inv.inner_control)
    #lcl
    p[19] = get_lf(inv.filter)
    p[20] = get_rf(inv.filter)
    p[21] = get_cf(inv.filter)
    p[22] = get_lg(inv.filter)
    p[23] = get_rg(inv.filter)
    return p
end

function build_train_system(sys_surr_original, sys_pvs_original, surrogate_area_name)
    sys_surr = deepcopy(sys_surr_original)
    sys_pvs = deepcopy(sys_pvs_original)
    non_surrogate_buses = collect(
        get_components(Bus, sys_surr, x -> get_name(get_area(x)) != surrogate_area_name),
    )
    if length(non_surrogate_buses) != 1
        @error "Must have one-non surrogate bus designated in surrogate system to add the PVS to"
        return
    end
    non_surrogate_bus = non_surrogate_buses[1]
    set_bustype!(non_surrogate_bus, BusTypes.REF)

    sources = get_components(Source, sys_pvs)
    for s in sources
        pvs = get_dynamic_injector(s)
        remove_component!(sys_pvs, pvs)
        remove_component!(sys_pvs, s)

        set_bus!(s, non_surrogate_bus)
        add_component!(sys_surr, s)
        add_component!(sys_surr, pvs, s)
    end
    return sys_surr
end

#Label an collection of buses with a name
function label_area!(sys::System, bus_numbers, area_name::String)
    buses = collect(get_components(Bus, sys))
    areas = collect(get_components(Area, sys))
    for area in areas
        if get_name(area) == area_name
            @error "area already exists"
            return 0
        end
    end
    surrogate_area = Area(; name = area_name)
    add_component!(sys, surrogate_area)
    for bus in buses
        if get_number(bus) in bus_numbers
            set_area!(bus, surrogate_area)
        end
    end
end

function check_single_connecting_line_condition(sys::System)
    areas = get_components(Area, sys)
    if length(areas) != 2
        @warn "There aren't two areas in this system. Cannot check single line condition"
    end
    branches = collect(get_components(Branch, sys))
    count_connections = 0
    for branch in branches
        area_from = get_name(get_area(get_from(get_arc(branch))))
        area_to = get_name(get_area(get_to(get_arc(branch))))
        if (area_from != area_to)
            count_connections += 1
        end
    end
    if count_connections == 1
        return true
    else
        @warn "Not exactly one connections between the two areas"
        return false
    end
end

function remove_area(sys_original::System, area_name::String)
    sys = deepcopy(sys_original)
    (length(collect(get_components(Area, sys, x -> get_name(x) == area_name))) == 1) ||
        @warn "area with name not found or multiple areas with same name"
    connecting_bus_name = nothing

    static_injectors =
        collect(get_components(Component, sys, x -> typeof(x) <: StaticInjection))
    for static_injector in static_injectors
        if get_name(get_area(get_bus(static_injector))) == area_name
            dynamic_injector = get_dynamic_injector(static_injector)
            (dynamic_injector !== nothing) && remove_component!(sys, dynamic_injector)
            remove_component!(sys, static_injector)
        end
    end

    branches = collect(get_components(Component, sys, x -> typeof(x) <: Branch))
    for branch in branches     #DOES ORDER MATTER ?  
        area_name_from = get_name(get_area(get_from(get_arc(branch))))
        area_name_to = get_name(get_area(get_to(get_arc(branch))))
        if (area_name_from == area_name_to == area_name)
            arc = get_arc(branch)
            remove_component!(sys, arc)
            remove_component!(sys, branch)
        end
        if (area_name_from != area_name_to)
            if (area_name_from == area_name)
                connecting_bus_name = get_name(get_from(get_arc(branch)))
            else
                connecting_bus_name = get_name(get_to(get_arc(branch)))
            end
        end
    end

    buses = collect(get_components(Component, sys, x -> typeof(x) <: Bus))
    for bus in buses
        if (get_name(get_area(bus)) == area_name) && (get_name(bus) != connecting_bus_name)
            remove_component!(sys, bus)
        end
    end
    return sys
end

function node_load_system(inputs...)
    logger =
        configure_logging(console_level = PSY_CONSOLE_LEVEL, file_level = PSY_FILE_LEVEL)
    try
        Logging.with_logger(logger) do
            sys = PSY.System(inputs...)
            return sys
        end
    finally
        close(logger)
        configure_logging(console_level = NODE_CONSOLE_LEVEL, file_level = NODE_FILE_LEVEL)
    end
end
