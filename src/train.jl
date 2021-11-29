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
