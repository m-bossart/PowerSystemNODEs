# TODO: Change name to _function to functions only used in this file for training.

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

function turn_node_on(surr_prob_node_off, params, p_ode, p_V₀, p_nn)
    n_weights_nn = length(p_nn)
    p_scale = [params.node_input_scale, params.node_output_scale]
    p_fixed = vcat(p_ode, p_scale, p_V₀, n_weights_nn)
    p = vcat(p_nn, p_fixed)
    surr_prob = remake(surr_prob_node_off, p = p)
    return surr_prob, p_fixed
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

function get_init_surr(p, ir_filter, ii_filter, surr)
    return (dx, x) -> begin
        dx[22] = 0
        x[22] = ir_filter
        dx[23] = 0
        x[23] = ii_filter
        surr(dx, x, p, 0.0)
    end
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

# TODO: We want to add types in here to make the function performant
function _train(
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
    sys = node_load_system(joinpath(params.input_data_path, "system.json"))
    TrainInputs =
        JSON3.read(read(joinpath(params.input_data_path, "data.json")), NODETrainInputs)
    d = TrainInputs.data
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

                res, output = _train(
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
