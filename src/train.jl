"""
# Fields
- `train_id = Int64`: id for the training instance, used for naming output data folder.
- `solver = ["Rodas4"]`: solver used for the NODE problem.
- `solver_tols =  Tuple{Float64, Float64}`: solver tolerances (abstol, reltol).
**Note:** Include a note 
"""
mutable struct NODETrainParams
    train_id::String
    solver::String
    solver_tols::Tuple{Float64, Float64}
    sensealg::String
    optimizer::String
    optimizer_η::Float64
    optimizer_adjust::String
    optimizer_adjust_η::Float64
    maxiters::Int64   #should this be for the full training or per call to the optimizer
    lb_loss::Float64  #default to 0 
    batching::Bool
    batch_factor::Float64 #multiply the size of the current training range by batch_factor to get the number of points.
    rng_seed::Int64     #re-create results?  #ONLY NEED ONE SEED 
    groupsize_steps::Int64  #the size of the extending timespan         
    groupsize_faults::Int64 #how many faults should be trained on simultaneously, default to 1. 
    loss_function_weights::Tuple{Float64, Float64}
    loss_function_scale::String     #"none", "range" 
    ode_model::String
    node_input_scale::Float64   #scale on the voltage input only
    node_output_scale::Float64  #scale on the current output only 
    node_inputs::String
    node_feedback_states::Int64
    node_feedback_current::Bool
    node_layers::Int64
    node_width::Int64
    node_activation::String
    output_mode::Int64
    base_path::String
    input_data_path::String
    output_data_path::String
    verify_psid_node_off::Bool
end

StructTypes.StructType(::Type{NODETrainParams}) = StructTypes.Struct()

function NODETrainParams(;
    train_id = "train_instance_1",
    solver = "Rodas4",
    solver_tols = (1e-6, 1e-9),
    sensealg = "ForwardDiffSensitivity",
    optimizer = "Adam",
    optimizer_η = 0.01,
    optimizer_adjust = "nothing",
    optimizer_adjust_η = 0.01,
    maxiters = 5,
    lb_loss = 0.0,
    batching = false,
    batch_factor = 1.0,
    rng_seed = 1234,
    groupsize_steps = 55,
    groupsize_faults = 1,
    loss_function_weights = (0.5, 0.5),
    loss_function_scale = "range",
    ode_model = "vsm",
    node_input_scale = 1e1,
    node_output_scale = 1.0,
    node_inputs = "voltage",
    node_feedback_states = 0,
    node_feedback_current = true,
    node_layers = 2,
    node_width = 2,
    node_activation = "relu",
    export_mode = 3,
    base_path = pwd(),
    input_data_path = joinpath(base_path, "input_data"),
    output_data_path = joinpath(base_path, "output_data"),
    verify_psid_node_off = true,
)

    #HERE IS THE LOGIC OF FILLING IN SOME OF THE PARAMETERS THAT MIGHT NOT MAKE SENSE       
    NODETrainParams(
        train_id,
        solver,
        solver_tols,
        sensealg,
        optimizer,
        optimizer_η,
        optimizer_adjust,
        optimizer_adjust_η,
        maxiters,
        lb_loss,
        batching,
        batch_factor,
        rng_seed,
        groupsize_steps,
        groupsize_faults,
        loss_function_weights,
        loss_function_scale,
        ode_model,
        node_input_scale,
        node_output_scale,
        node_inputs,
        node_feedback_states,
        node_feedback_current,
        node_layers,
        node_width,
        node_activation,
        export_mode,
        base_path,
        input_data_path,
        output_data_path,
        verify_psid_node_off,
    )
end

function read_input_data(pvs, d)
    id = get_name(pvs)
    tsteps = Float64.(d[id][:tsteps])
    i_ground_truth =
        vcat(Float64.(d[id][:ir_ground_truth])', Float64.(d[id][:ii_ground_truth])')
    i_node_off = vcat(Float64.(d[id][:ir_node_off])', Float64.(d[id][:ii_node_off])')
    p_ode = Float64.(d[id][:p_ode])
    x₀ = Float64.(d[id][:x₀])
    p_V₀ = Float64.(d[id][:V₀])
    return id, tsteps, i_ground_truth, i_node_off, p_ode, x₀, p_V₀
end

function train(params::NODETrainParams)
    #INSTANTIATE 
    sensealg = instantiate_sensealg(params)
    solver = instantiate_solver(params)
    optimizer = instantiate_optimizer(params)
    nn = instantiate_nn(params)
    M = instantiate_M(params)
    !(params.optimizer_adjust == "nothing") &&
        (optimizer_adjust = instantiate_optimizer_adjust(params))

    #READ INPUT DATA AND SYSTEM
    sys = System(joinpath(params.input_data_path, "system.json"))
    d = JSON3.read(
        read(joinpath(params.input_data_path, "data.json")),
        Dict{String, Dict{Symbol, Any}},
    )
    pvss = collect(get_components(PeriodicVariableSource, sys))

    #TRAIN SEQUENTIALLY 
    res = nothing
    output = Dict{String, Any}(
        "loss" => DataFrame(ID = String[], RangeCount = Int[], Loss = Float64[]),
        "parameters" => DataFrame(Parameters = Vector{Any}[]),
        "predictions" =>
            DataFrame(ir_prediction = Vector{Any}[], ii_prediction = Vector{Any}[]),
        "total_time" => [],
        "total_iterations" => [],
        "final_results" => Dict{String, Any}(),
    )

    min_θ = initial_params(nn)
    for pvs in pvss
        @show min_θ[end]
        res, output =
            train(min_θ, params, sensealg, solver, optimizer, nn, M, d, pvs, output)
        min_θ = copy(res.u)
    end

    #TRAIN ADJUSTMENTS (TO DO)
    @show output
    capture_output(output, params.output_data_path, params.train_id)

    return output
end

function train(θ, params, sensealg, solver, optimizer, nn, M, d, pvs, output) #move 159-168 outside of train
    #READ FAULT DATA FOR THE CURRENT PVS 
    id, tsteps, i_true, i_ver, p_ode, x₀, p_V₀ = read_input_data(pvs, d)

    #DEFINE FUNCTIONS OF TIME FOR THE CURRENT PVS
    Vm, Vθ = Source_to_function_of_time(pvs)
    surr = instantiate_surr(params, nn, Vm, Vθ)

    #DEFINE SCALE FACTORS FOR THE CURRENT PVS , TODO  REMOVE 
    Vr_scale =
        1 / (
            maximum(Vm.(tsteps) .* cos.(Vθ.(tsteps))) -
            minimum(Vm.(tsteps) .* cos.(Vθ.(tsteps)))
        )
    Vi_scale =
        1 / (
            maximum(Vm.(tsteps) .* sin.(Vθ.(tsteps))) -
            minimum(Vm.(tsteps) .* sin.(Vθ.(tsteps)))
        )

    #FIND INITIAL CONDITIONS FOR THE SURROGATE, BUILD AND SOLVE THE PROBLEM, AND CONFIRM IT MATCHES THE VERIFICATION DATA PROVIDED. ABSTRACT TO FUNCITON (170-191), call inside train 
    p_scale = [Vr_scale, Vi_scale, 0.0]     #turn off the nn 
    p_nn = initial_params(nn)
    n_weights_nn = length(p_nn)
    p_fixed = vcat(p_ode, p_scale, p_V₀, n_weights_nn)
    p = vcat(p_nn, p_fixed)
    order_surr = ODE_ORDER + 2 + params.node_feedback_states + 2   #2 states for the node current, 2 algebraic states   
    x₀_surr = zeros(order_surr)
    x₀_surr[1:ODE_ORDER] = x₀
    x₀_surr[(end - 1):end] = [x₀[I__IR_FILTER], x₀[I__II_FILTER]]
    h = get_init_vsm_v_t_0(p, x₀[I__IR_FILTER], x₀[I__II_FILTER], nn, Vm, Vθ) #make generic! get rid of hard code 
    res_surr = nlsolve(h, x₀_surr)
    @assert converged(res_surr)
    dx = similar(x₀_surr)
    surr(dx, res_surr.zero, p, 0.0)
    @assert all(isapprox.(dx, 0.0; atol = 1e-8))

    tspan = (tsteps[1], tsteps[end])   #tspan comes from data   #Pass flag as parameter in NODETrainParams if you want to verify or not 
    surr_func = ODEFunction(surr, mass_matrix = M)
    surr_prob = ODEProblem(surr_func, x₀_surr, tspan, p)
    sol = solve(
        surr_prob,
        solver,
        abstol = params.solver_tols[1],
        reltol = params.solver_tols[2],
        saveat = tsteps,
    )
    @assert mae(sol[22, :], i_ver[1, :]) < 5e-5

    #PREPARE THE SURROGATE FOR TRAINING  
    p_scale = [Vr_scale, Vi_scale, params.node_output_scale]
    p_fixed = vcat(p_ode, p_scale, p_V₀, n_weights_nn)
    p = vcat(p_nn, p_fixed)
    surr_prob = remake(surr_prob, p = p)

    #Calculate scale factors for loss function based on true data - Create a function for all of the scaling. 
    Ir_scale = maximum(i_true[1, :]) - minimum(i_true[1, :])
    Ii_scale = maximum(i_true[2, :]) - minimum(i_true[2, :])
    #@show surr_init = instantiate_surr_init(params) #don't need? just use surr? 

    u₀ = res_surr.zero

    #INSTANTIATE THE TRAINING FUNCTIONS 
    pred_function = instantiate_pred_function(
        p_fixed,
        solver,
        surr_prob,
        params.solver_tols,
        sensealg,
        u₀,
    )
    loss_function = instantiate_loss_function(
        params.loss_function_weights,
        Ir_scale,
        Ii_scale,
        pred_function,
    )

    #TRAIN ON A SINGLE FAULT USING EXTENDING TIME RANGES.
    datasize = length(tsteps)
    ranges = extending_ranges(datasize, params.groupsize_steps)
    res = nothing
    min_θ = θ
    range_count = 1
    for range in ranges
        @show min_θ[end]
        i_curr = i_true[:, range]
        t_curr = tsteps[range]
        train_loader = Flux.Data.DataLoader(
            (i_curr, t_curr),
            batchsize = Int(floor(length(i_curr[1, :]))),
        )     #TODO - IMPLEMENT BATCHING
        optfun = OptimizationFunction(
            (θ, p, batch, time_batch) -> loss_function(θ, batch, time_batch),
            GalacticOptim.AutoForwardDiff(),
        )
        optprob = OptimizationProblem(optfun, min_θ)

        cb = instantiate_cb!(output, params.lb_loss, params.output_mode, id, range_count)
        range_count += 1

        res = GalacticOptim.solve(
            optprob,
            optimizer,
            ncycle(train_loader, params.maxiters),
            cb = cb,
        )
        min_θ = copy(res.u)
        @show min_θ[end]
        @assert res.minimum == loss_function(res.u, i_curr, t_curr)[1]
        @assert res.minimum == loss_function(min_θ, i_curr, t_curr)[1]
    end
    @show min_θ
    return res, output
    try
    catch e
        return 0
    end
end

function capture_output(output_dict, output_directory, id)
    output_path = joinpath(output_directory, id)
    mkpath(output_path)
    for (key, value) in output_dict
        if typeof(value) == DataFrame
            df = pop!(output_dict, key)
            open(joinpath(output_path, key), "w") do io
                Arrow.write(io, df)
            end
        end
    end
    open(joinpath(output_path, "high_level_outputs"), "w") do io
        JSON3.write(io, output_dict)
    end
end
