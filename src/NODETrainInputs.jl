# TODO: Use Dict temporarily during dev while the fields are defined
struct NODETrainInputs
    inputs_name::String
    data::Dict{Symbol, Vector{Float64}}
end

function NODETrainInputs(name::String)
    return NODETrainInputs(name, Dict{Symbol, Vector{Float64}}())
end

function serialize(inputs::NODETrainInputs, file_path::String)
    open(file_path, "w") do io
        JSON3.write(io, inputs)
    end
    return
end

mutable struct NODETrainDataParams
    solver::String
    solver_tols::Tuple{Float64, Float64}
    tspan::Tuple{Float64, Float64}
    steps::Int64
    tsteps_spacing::String
    base_path::String
    output_data_path::String
end

StructTypes.StructType(::Type{NODETrainDataParams}) = StructTypes.Struct()
StructTypes.StructType(::Type{NODETrainInputs}) = StructTypes.Struct()
function NODETrainDataParams(;
    solver = "Rodas4",
    solver_tols = (1e-6, 1e-9),
    tspan = (0.0, 2.0),
    steps = 150,
    tsteps_spacing = "linear",
    base_path = pwd(),
    output_data_path = joinpath(base_path, "input_data"),
)
    NODETrainDataParams(
        solver,
        solver_tols,
        tspan,
        steps,
        tsteps_spacing,
        base_path,
        output_data_path,
    )
end

function generate_train_data(sys_train, NODETrainDataParams)
    tspan = NODETrainDataParams.tspan
    steps = NODETrainDataParams.steps
    if NODETrainDataParams.tsteps_spacing == "linear"
        tsteps = tspan[1]:((tspan[2] - tspan[1]) / steps):tspan[2]
    end
    solver = instantiate_solver(NODETrainDataParams)
    abstol = NODETrainDataParams.solver_tols[1]
    reltol = NODETrainDataParams.solver_tols[2]

    available_source = activate_next_source!(sys_train)

    set_bus_from_source(available_source) #Bus voltage is used in power flow, not source voltage. Need to set bus voltage from soure internal voltage

    sim_full = Simulation!(MassMatrixModel, sys_train, pwd(), tspan)
    #res = small_signal_analysis(sim_full)
    execute!(
        sim_full,
        solver,
        abstol = abstol,
        reltol = reltol,
        reset_simulation = false,
        saveat = tsteps,
    )
    active_source = collect(get_components(Source, sys_train, x -> PSY.get_available(x)))[1]
    ode_data = get_total_current_series(sim_full) #TODO Better to measure current at the PVS (implement method after PVS is complete)

    #################### BUILD INITIALIZATION SYSTEM ###############################
    sys_init, p_inv = build_sys_init(sys_train) #returns p_inv, the set of average parameters 
    transformer = collect(get_components(Transformer2W, sys_init))[1]
    pvs = collect(get_components(PeriodicVariableSource, sys_init))[1]
    p_fixed = [get_x(transformer) + get_X_th(pvs), get_r(transformer) + get_R_th(pvs)]
    x₀, refs, Vr0, Vi0 = initialize_sys!(sys_init, "gen1")
    Vm, Vθ = Source_to_function_of_time(get_dynamic_injector(active_source))
    p_ode = vcat(p_inv, refs, p_fixed)
    sim_simp = Simulation!(MassMatrixModel, sys_init, pwd(), tspan)
    @debug "initialize system power flow", solve_powerflow(sys_init)["flow_results"]
    @debug "initialize system power flow", solve_powerflow(sys_init)["bus_results"]
    @debug show_states_initial_value(sim_simp)
    @time execute!(
        sim_simp,
        solver,
        abstol = abstol,
        reltol = reltol,
        initializealg = NoInit(),
        reset_simulation = false,
        saveat = tsteps,
    )

    #avgmodel_data_p = get_real_current_series(read_results(sim_simp), "gen1")
    avgmodel_data = get_total_current_series(sim_simp)

    d = NODETrainInputs(
        get_name(pvs),
        Dict(
            :tsteps => tsteps,
            :ir_ground_truth => ode_data[1, :],
            :ii_ground_truth => ode_data[2, :],
            :ir_node_off => avgmodel_data[1, :],
            :ii_node_off => avgmodel_data[2, :],
            :p_ode => p_ode,
            :x₀ => x₀,
            :V₀ => [Vr0, Vi0],
        ),
    )

    return d
end
