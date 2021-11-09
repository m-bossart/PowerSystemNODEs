"""
# Fields
- `train_id = Int64`: id for the training instance, used for naming output data folder.
- `solver = ["Rodas4"]`: solver used for the NODE problem.
- `solver_tols =  Tuple{Float64, Float64}`: solver tolerances (abstol, reltol).
- `sensealg = ["ForwardDiffSensitivity"]`: sensitivity algorithm used in training. 
- `optimizer = ["Adam", "Bfgs"]`: main optimizer used in training. 
- `optimizer_η = Float64`: Learning rate for Adam (amount by which gradients are discounted before updating weights). Ignored if Adam is not the optimizer.
- `optimizer_adjust = ["Adam", "Bfgs", "nothing"]`: optimizer used for final adjustments (2nd stage). 
- `optimizer_adjust_η = Float64`: Learning rate for Adam (amount by which gradients are discounted before updating weights). Ignored if Adam is not the optimizer.
- `maxiters = Int64`: The maximum possible iterations for the entire training instance. If `lb_loss = 0` and `optimizer = "Adam"` the training should never exit early and maxiters will be hit. 
    Note that the number of saved data points can exceed maxiters because there is an additional callback at the end of each individual optimization.  
- `lb_loss = Float64`: If the value of the loss function moves below lb_loss during training, the current optimization ends (current range).
- `batching = Bool`: If `batching = false` the full set of data points are used for each training step.
- `batching_factor = Float64`: The number of data points in the current range is multiplied by `batching_factor` to get the size of the batch. Batches of this size are used sequentially in time. 
    The final batch is used even if it is incomplete.  
**Note:** BATCHING IS NOT YET IMPLEMENTED
- `rng_seed = Int64`: Seed for the random number generator used for initializing the NN for reproducibility across training runs.  
- `groupsize_steps = Int64`: Number of data-points in each extension of the range of data used. 
- `groupsize_faults = Int64`: Number of data-points in each extension of the range of data used. 
**Note:** GROUPSIZE_FAULTS NOT YET IMPLEMENTED. NOT NEEDED FOR SINGLE FAULT TRAINING.
- `loss_function_weights = Tuple{Float64, Float64}`: weights used for loss function `(mae_weight, mse_weight)`. 
- `loss_function_scale = ["range", "none"]`: Scaling of the loss function.  `"range"`: the range of the real current and imaginary current are used to scale both the mae  
    and mse portions of the loss function. The goal is to give equal weight to real and imaginary components even if the magnitude of the disturbance differs. `"none"`: no additional scaling applied.
- `ode_model = ["vsm"]`: The ode model used in conjunction with the NODE during training.
- `node_input_scale = Float64`: Scale factor on the voltage input to the NODE. Does not apply to other inputs (ie the feedback states).
- `node_output_scale = Float64`: Scale factor on the current output of the NODE. Does not apply to other outputs (ie the feedback states).
- `node_inputs = ["voltage"]`: Determines the physical states which are inputs to the NODE. Ideally, only voltage to remain as general as possible. 
- `node_feedback_states = Int64`: Number of feedback states in the NODE. Does not include the output current states which can be feedback if `node_feedback_current = true`. 
- `node_feedback_current = Bool`: Determines if current is also a feedback state. 
- `node_layers = Int64`: Number of hidden layers in the NODE. Does not include the input or output layer.
- `node_width = Int64`: Number of neurons in each hidden layer. Each hidden layer has the same number of neurons. The width of the input and output layers are determined by the combination of other parameters. 
- `node_activation = ["relu"]`: Activation function for NODE. The output layer always uses the identity activation.  
- `output_mode = [1,2,3]`: `1`: do not collect any data during training, only save high-level data related to training and final results `2`: Same as `1`, also save value of loss throughout training.
    `3`: same as `2`, also save parameters and predictions during training. 
- `base_path = String`: Directory for training where input data is found and output data is written. 
- `input_data_path = String`: From `base_path`, the directory for input data.
- `output_data_path = String`: From `base_path`, the directory for saving output data.
- `verify_psid_node_off = Bool`: `true`: before training, check that the surrogate with NODE turned off matches the data provided from PSID simulation. 
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
    maxiters::Int64
    lb_loss::Float64
    batching::Bool
    batch_factor::Float64
    rng_seed::Int64
    groupsize_steps::Int64
    groupsize_faults::Int64
    loss_function_weights::Tuple{Float64, Float64}
    loss_function_scale::String
    ode_model::String
    node_input_scale::Float64
    node_output_scale::Float64
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
    maxiters = 15,
    lb_loss = 0.0,
    batching = false,
    batch_factor = 1.0,
    rng_seed = 1234,
    groupsize_steps = 55,
    groupsize_faults = 1,
    loss_function_weights = (0.5, 0.5),
    loss_function_scale = "range",
    ode_model = "vsm",
    node_input_scale = 10e1,
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
function calculate_per_solve_maxiters(params, pvss, d)
    n_faults = length(pvss)
    id, tsteps, i_true, i_ver, p_ode, x₀, p_V₀ = read_input_data(pvss[1], d)
    n_timesteps = length(tsteps)
    total_maxiters = params.maxiters
    groupsize_steps = params.groupsize_steps
    factor_ranges = ceil(n_timesteps / groupsize_steps)
    factor_batches = ceil(1 / params.batch_factor)
    per_solve_maxiters =
        Int(floor(total_maxiters / factor_ranges / factor_batches / n_faults))
    @info "per solve maxiters" per_solve_maxiters
    if per_solve_maxiters == 0
        @error "maxiters is too low. The calculated maxiters per solve is 0! cannot train"
    end
    return per_solve_maxiters
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

    res = nothing
    output = Dict{String, Any}(
        "loss" => DataFrame(ID = String[], RangeCount = Int[], Loss = Float64[]),
        "parameters" => DataFrame(Parameters = Vector{Any}[]),
        "predictions" =>
            DataFrame(ir_prediction = Vector{Any}[], ii_prediction = Vector{Any}[]),
        "total_time" => [],
        "total_iterations" => 0,
        "final_loss" => [],
        "train_id" => params.train_id,
    )
    per_solve_maxiters = calculate_per_solve_maxiters(params, pvss, d)
    min_θ = initial_params(nn)
    try
        total_time = @elapsed begin
            #TRAIN SEQUENTIALLY 
            for pvs in pvss
                @info "start of fault" min_θ[end]

                id, tsteps, i_true, i_ver, p_ode, x₀, p_V₀ = read_input_data(pvs, d)

                Vm, Vθ = Source_to_function_of_time(pvs)
                surr = instantiate_surr(params, nn, Vm, Vθ)

                res, output = train(
                    min_θ,
                    params,
                    sensealg,
                    solver,
                    optimizer,
                    nn,
                    M,
                    output,
                    id,
                    tsteps,
                    i_true,
                    i_ver,
                    p_ode,
                    x₀,
                    p_V₀,
                    surr,
                    Vm,
                    Vθ,
                    per_solve_maxiters,
                )
                min_θ = copy(res.u)
                @info "end of fault" min_θ[end]
            end

            #TRAIN ADJUSTMENTS GO HERE (TO DO)
        end
        @info "min_θ[end] (end of training)" min_θ[end]
        output["total_time"] = total_time

        final_loss_for_comparison =
            calculate_final_loss(params, res.u, solver, nn, M, pvss, d, sensealg)
        output["final_loss"] = final_loss_for_comparison

        capture_output(output, params.output_data_path, params.train_id)
        return true
    catch
        return false
    end
end

function calculate_final_loss(params, θ, solver, nn, M, pvss, d, sensealg)
    final_loss_for_comparison = 0.0
    n_pvs = length(pvss)
    for pvs in pvss
        id, tsteps, i_true, i_ver, p_ode, x₀, p_V₀ = read_input_data(pvs, d)
        Vm, Vθ = Source_to_function_of_time(pvs)
        surr = instantiate_surr(params, nn, Vm, Vθ)

        u₀, surr_prob_node_off, p_nn, p_fixed =
            initialize_surrogate(params, nn, M, tsteps, p_ode, x₀, p_V₀, surr, Vm, Vθ)

        surr_prob, p_fixed = turn_node_on(surr_prob_node_off, params, p_ode, p_V₀, p_nn)

        pred_function = instantiate_pred_function(
            p_fixed,
            solver,
            surr_prob,
            params.solver_tols,
            sensealg,
            u₀,
        )
        loss_function = instantiate_loss_function(
            (1.0, 0.0), #mae only
            1.0,        #no scaling
            1.0,
            pred_function,
        )
        final_loss_for_comparison += loss_function(θ, i_true, tsteps)[1]
    end

    return final_loss_for_comparison / n_pvs
end

function verify_psid_node_off(surr_prob, params, solver, tsteps, i_ver)
    sol = solve(
        surr_prob,
        solver,
        abstol = params.solver_tols[1],
        reltol = params.solver_tols[2],
        saveat = tsteps,
    )
    @assert mae(sol[22, :], i_ver[1, :]) < 5e-5
end

function initialize_surrogate(params, nn, M, tsteps, p_ode, x₀, p_V₀, surr, Vm, Vθ)
    p_scale = [params.node_input_scale, 0.0]     #turn off the nn 
    p_nn = initial_params(nn)
    n_weights_nn = length(p_nn)
    p_fixed = vcat(p_ode, p_scale, p_V₀, n_weights_nn)
    p = vcat(p_nn, p_fixed)
    order_surr = ODE_ORDER + 2 + params.node_feedback_states + 2   #2 states for the node current, 2 algebraic states   
    x₀_surr = zeros(order_surr)
    x₀_surr[1:ODE_ORDER] = x₀
    x₀_surr[(end - 1):end] = [x₀[I__IR_FILTER], x₀[I__II_FILTER]]
    h = get_init_surr(p, x₀[I__IR_FILTER], x₀[I__II_FILTER], surr) #make generic! get rid of hard code 
    res_surr = nlsolve(h, x₀_surr)
    @assert converged(res_surr)
    dx = similar(x₀_surr)
    surr(dx, res_surr.zero, p, 0.0)
    @assert all(isapprox.(dx, 0.0; atol = 1e-8))

    tspan = (tsteps[1], tsteps[end])   #tspan comes from data   #Pass flag as parameter in NODETrainParams if you want to verify or not 
    surr_func = ODEFunction(surr, mass_matrix = M)
    surr_prob = ODEProblem(surr_func, x₀_surr, tspan, p)
    return res_surr.zero, surr_prob, p_nn, p_fixed
end

function turn_node_on(surr_prob_node_off, params, p_ode, p_V₀, p_nn)
    n_weights_nn = length(p_nn)
    p_scale = [params.node_input_scale, params.node_output_scale]
    p_fixed = vcat(p_ode, p_scale, p_V₀, n_weights_nn)
    p = vcat(p_nn, p_fixed)
    surr_prob = remake(surr_prob_node_off, p = p)
    return surr_prob, p_fixed
end

function calculate_loss_function_scaling(params, i_true)
    if params.loss_function_scale == "range"
        Ir_scale = maximum(i_true[1, :]) - minimum(i_true[1, :])
        Ii_scale = maximum(i_true[2, :]) - minimum(i_true[2, :])
    elseif params.loss_function_scale == "none"
        Ir_scale = 1.0
        Ii_scale = 1.0
    else
        @warn "Cannot determine loss function scaling"
    end
    return Ir_scale, Ii_scale
end

function train(
    θ,
    params,
    sensealg,
    solver,
    optimizer,
    nn,
    M,
    output,
    id,
    tsteps,
    i_true,
    i_ver,
    p_ode,
    x₀,
    p_V₀,
    surr,
    Vm,
    Vθ,
    per_solve_maxiters,
)
    u₀, surr_prob_node_off, p_nn, p_fixed =
        initialize_surrogate(params, nn, M, tsteps, p_ode, x₀, p_V₀, surr, Vm, Vθ)

    (params.verify_psid_node_off) &&
        verify_psid_node_off(surr_prob_node_off, params, solver, tsteps, i_ver)

    surr_prob, p_fixed = turn_node_on(surr_prob_node_off, params, p_ode, p_V₀, p_nn)

    Ir_scale, Ii_scale = calculate_loss_function_scaling(params, i_true)

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
        @info "start of range" min_θ[end]
        i_curr = i_true[:, range]
        t_curr = tsteps[range]
        batchsize = Int(floor(length(i_curr[1, :]) * params.batch_factor))
        train_loader = Flux.Data.DataLoader(
            (i_curr, t_curr),
            batchsize = batchsize,   #TODO - IMPLEMENT BATCHING
        )
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
            ncycle(train_loader, per_solve_maxiters),
            cb = cb,
        )
        min_θ = copy(res.u)
        @info "end of range" min_θ[end]
        if params.batch_factor == 1.0
            @assert res.minimum == loss_function(res.u, i_curr, t_curr)[1]
            @assert res.minimum == loss_function(min_θ, i_curr, t_curr)[1]
        end
    end
    return res, output
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

function get_init_surr(p, ir_filter, ii_filter, surr)
    return (dx, x) -> begin
        dx[22] = 0
        x[22] = ir_filter
        dx[23] = 0
        x[23] = ii_filter
        surr(dx, x, p, 0.0)
    end
end
