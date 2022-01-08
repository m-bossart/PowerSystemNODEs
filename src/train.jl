# TODO: Change name to _function to functions only used in this file for training.

function calculate_loss_function_scaling(params, fault_data)
    if params.loss_function_scale == "range"
        full_ir = Float64[]
        full_ii = Float64[]
        for (key, value) in fault_data
            full_ir = vcat(full_ir, value[:ir_ground_truth])
            full_ii = vcat(full_ii, value[:ii_ground_truth])
        end
        Ir_scale = maximum(full_ir) - minimum(full_ir)
        Ii_scale = maximum(full_ii) - minimum(full_ii)
    elseif params.loss_function_scale == "none"
        Ir_scale = 1.0
        Ii_scale = 1.0
    else
        @warn "Cannot determine loss function scaling"
    end
    return Ir_scale, Ii_scale
end
function turn_node_on(surr_prob_node_off, params, P)
    P.scale[2] = params.node_output_scale
    p = vectorize(P)
    surr_prob = remake(surr_prob_node_off, p = p)
    return surr_prob
end

function initialize_surrogate(
    params,
    nn,
    M,
    tsteps,
    fault_dict,
    surr,
    N_ALGEBRAIC_STATES,
    ODE_ORDER,
)
    #Determine Order of the Surrogate 
    order_surr = ODE_ORDER + 2 + params.node_feedback_states + N_ALGEBRAIC_STATES    #2 states for the node current
    x₀_surr = zeros(order_surr)

    #Build Surrogate Vector 
    P = SurrParams()
    if params.ode_model == "none"
        if !(isempty(fault_dict[:p_ode]))
            @warn "ODE model is not present in surrogate, yet parameters for an ode are provided in train data. These parameters will be ignored"
        end
        P.ode = []
    else
        P.ode = fault_dict[:p_ode]
    end
    P.pf = fault_dict[:p_pf]
    P.network = fault_dict[:p_network]
    P.nn = initial_params(nn)
    P.n_weights = [Float64(length(P.nn))]
    if params.ode_model == "none"
        P.scale = [params.node_input_scale, params.node_output_scale]
        #x₀_surr[1:(end-2)] = 0.0  - initialized to zero already 
    else
        P.scale = [params.node_input_scale, 0.0]
        x₀ = fault_dict[:x₀]
        x₀_surr[(3 + params.node_feedback_states + N_ALGEBRAIC_STATES):(ODE_ORDER + 2 + params.node_feedback_states + N_ALGEBRAIC_STATES)] =
            x₀
    end
    p = vectorize(P)

    P_pf, Q_pf, V_pf, θ_pf = P.pf
    Vr_pf = V_pf * cos(θ_pf)
    Vi_pf = V_pf * sin(θ_pf)
    Ir_pf = (P_pf * Vr_pf + Q_pf * Vi_pf) / (Vr_pf^2 + Vi_pf^2)
    Ii_pf = (P_pf * Vi_pf - Q_pf * Vr_pf) / (Vi_pf^2 + Vr_pf^2)

    x₀_surr[1:2] = [Ir_pf, Ii_pf]

    if params.ode_model != "none"   #only do an NL solve if there is an ODE part. Otherwise, set initial conditions and try to solve.
        @warn "Power flow results", P.pf
        h = get_init_surr(p, Ir_pf, Ii_pf, surr)
        res_surr = nlsolve(h, x₀_surr)
        @assert converged(res_surr)
        dx = similar(x₀_surr)
        surr(dx, res_surr.zero, p, 0.0)
        @assert all(isapprox.(dx, 0.0; atol = 1e-8))
        x₀_surr = res_surr.zero
    end
    @warn "x0 surr", x₀_surr
    tspan = (tsteps[1], tsteps[end])
    surr_func = ODEFunction(surr, mass_matrix = M)
    surr_prob = ODEProblem(surr_func, x₀_surr, tspan, p)
    return surr_prob, P
end

function verify_psid_node_off(surr_prob, params, solver, tsteps, fault_dict)
    i_ver = vcat(
        Float64.(fault_dict[:ir_ground_truth])',
        Float64.(fault_dict[:ii_ground_truth])',
    )
    sol = solve(
        surr_prob,
        solver,
        abstol = params.solver_tols[1],
        reltol = params.solver_tols[2],
        saveat = tsteps,
    )
    @show mae(sol[22, :], i_ver[1, :])
    @assert mae(sol[22, :], i_ver[1, :]) < 1e-3 # was 5e-5 with sequential train. need to double check visually 
end

function calculate_final_loss(
    params,
    θ,
    solver,
    nn,
    M,
    pvs_names,
    fault_data,
    tsteps,
    sensealg,
    Ir_scale,
    Ii_scale,
    output,
)
    t_concatonated = concatonate_t(tsteps, pvs_names, :)

    cb = instantiate_cb!(
        output,
        params.lb_loss,
        params.output_mode,
        0,  #range_count
        pvs_names,
        t_concatonated,
    )
    pred_function = instantiate_pred_function(
        solver,
        pvs_names,
        fault_data,
        params.solver_tols,
        sensealg,
    )
    loss_function = instantiate_loss_function(
        params.loss_function_weights,
        Ir_scale,
        Ii_scale,
        pred_function,
    )
    i_true = concatonate_i_true(fault_data, pvs_names, :)
    t_current = concatonate_t(tsteps, pvs_names, :)
    pvs_names = concatonate_pvs_names(pvs_names, length(tsteps))
    final_loss_for_comparison = loss_function(θ, i_true, t_current, pvs_names)

    cb(θ, final_loss_for_comparison[1], final_loss_for_comparison[2])     #Call the callback to record the prediction in output
    return final_loss_for_comparison[1]
end

function get_init_surr(p, Ir_pf, Ii_pf, surr)
    return (dx, x) -> begin
        dx[1] = 0
        x[1] = Ir_pf
        dx[2] = 0
        x[2] = Ii_pf
        surr(dx, x, p, 0.0)
    end
end

function calculate_per_solve_maxiters(params, tsteps, n_faults)
    n_timesteps = length(tsteps)
    total_maxiters = params.maxiters
    groupsize_steps = params.groupsize_steps
    groupsize_faults = params.groupsize_faults
    factor_ranges = ceil(n_timesteps / groupsize_steps)
    factor_faults = ceil(groupsize_faults)
    @warn factor_faults
    factor_batches = ceil(1 / params.batch_factor)
    per_solve_maxiters = Int(
        floor(total_maxiters * factor_faults / factor_ranges / factor_batches / n_faults),
    )
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
    sys = node_load_system(joinpath(params.input_data_path, "system.json"))

    TrainInputs =
        JSON3.read(read(joinpath(params.input_data_path, "data.json")), NODETrainInputs)

    tsteps = TrainInputs.tsteps
    fault_data = TrainInputs.fault_data
    pvss = collect(get_components(PeriodicVariableSource, sys))
    Ir_scale, Ii_scale = calculate_loss_function_scaling(params, fault_data)

    res = nothing
    output = Dict{String, Any}(
        "loss" => DataFrame(
            PVS_name = Vector{String}[],
            RangeCount = Int[],
            Loss = Float64[],
        ),
        "parameters" => DataFrame(Parameters = Vector{Any}[]),
        "predictions" => DataFrame(
            t_prediction = Vector{Any}[],
            ir_prediction = Vector{Any}[],
            ii_prediction = Vector{Any}[],
        ),
        "total_time" => [],
        "total_iterations" => 0,
        "final_loss" => [],
        "train_id" => params.train_id,
    )
    per_solve_maxiters =
        calculate_per_solve_maxiters(params, TrainInputs.tsteps, length(pvss))

    #PREPARE SURROGATES FOR EACH FAULT - this is what is different for pure NODE 
    for pvs in pvss
        Vm, Vθ = Source_to_function_of_time(pvs)
        surr, N_ALGEBRAIC_STATES, ODE_ORDER = instantiate_surr(params, nn, Vm, Vθ)
        fault_dict = fault_data[get_name(pvs)]
        surr_prob_node_off, P = initialize_surrogate(   #change naming
            params,
            nn,
            M,
            tsteps,
            fault_dict,
            surr,
            N_ALGEBRAIC_STATES,
            ODE_ORDER,
        )

        if (params.ode_model != "none")
            (params.verify_psid_node_off) &&
                verify_psid_node_off(surr_prob_node_off, params, solver, tsteps, fault_dict)
            surr_prob = turn_node_on(surr_prob_node_off, params, P)
        else
            surr_prob = surr_prob_node_off
        end

        fault_data[get_name(pvs)][:surr_problem] = surr_prob
        fault_data[get_name(pvs)][:P] = P   #parameters are stored in surr_prob, but as a single vector. The P struct allows for easier handling. 
    end

    min_θ = initial_params(nn)
    try
        total_time = @elapsed begin
            for group_pvs in partition(pvss, params.groupsize_faults)
                @info "start of fault" min_θ[end]
                @show pvs_names_subset = get_name.(group_pvs)

                res, output = _train(
                    min_θ,
                    params,
                    sensealg,
                    solver,
                    optimizer,
                    Ir_scale,
                    Ii_scale,
                    output,
                    tsteps,
                    pvs_names_subset,
                    fault_data,
                    per_solve_maxiters,
                )

                min_θ = copy(res.u)
                @info "end of fault" min_θ[end]
            end

            #TRAIN ADJUSTMENTS GO HERE (TODO)
        end
        @info "min_θ[end] (end of training)" min_θ[end]
        output["total_time"] = total_time

        pvs_names = get_name.(pvss)
        final_loss_for_comparison = calculate_final_loss(
            params,
            res.u,
            solver,
            nn,
            M,
            pvs_names,
            fault_data,
            tsteps,
            sensealg,
            Ir_scale,
            Ii_scale,
            output,
        )
        output["final_loss"] = final_loss_for_comparison

        capture_output(output, params.output_data_path, params.train_id)
        (params.graphical_report_mode != 0) &&
            visualize_training(params, visualize_level = params.graphical_report_mode)
        return true
    catch
        return false
    end
end

# TODO: We want to add types in here to make the function performant
function _train(
    θ,
    params,
    sensealg,
    solver,
    optimizer,
    Ir_scale,
    Ii_scale,
    output,
    tsteps,
    pvs_names_subset,
    fault_data,
    per_solve_maxiters,
)
    pred_function = instantiate_pred_function(
        solver,
        pvs_names_subset,
        fault_data,
        params.solver_tols,
        sensealg,
    )

    loss_function = instantiate_loss_function(
        params.loss_function_weights,
        Ir_scale,
        Ii_scale,
        pred_function,
    )

    datasize = length(tsteps)
    ranges = extending_ranges(datasize, params.groupsize_steps)
    res = nothing
    min_θ = θ
    range_count = 1
    for range in ranges
        @info "start of range" min_θ[end]
        i_current_range = concatonate_i_true(fault_data, pvs_names_subset, range)
        t_current_range = concatonate_t(tsteps, pvs_names_subset, range)
        pvs_names_current_range = concatonate_pvs_names(pvs_names_subset, length(range))
        batchsize = Int(floor(length(i_current_range[1, :]) * params.batch_factor))
        train_loader = Flux.Data.DataLoader(
            (i_current_range, t_current_range, pvs_names_current_range),
            batchsize = batchsize,
        )   #TODO - default for shuffle is false, make new parameter? 

        optfun = OptimizationFunction(
            (θ, p, batch, time_batch, pvs_name_batch) ->
                loss_function(θ, batch, time_batch, pvs_name_batch),
            GalacticOptim.AutoForwardDiff(),
        )
        optprob = OptimizationProblem(optfun, min_θ)
        cb = instantiate_cb!(
            output,
            params.lb_loss,
            params.output_mode,
            range_count,
            pvs_names_subset,
            t_current_range,
        )
        range_count += 1

        res = GalacticOptim.solve(
            optprob,
            optimizer,
            ncycle(train_loader, per_solve_maxiters),
            cb = cb,
        )
        min_θ = copy(res.u)
        @info "end of range" min_θ[end]
        if params.batch_factor == 1.0   #check that the minimum value is returned as expected. 
            @assert res.minimum == loss_function(
                res.u,
                i_current_range,
                t_current_range,
                pvs_names_current_range,
            )[1]
            @assert res.minimum == loss_function(
                min_θ,
                i_current_range,
                t_current_range,
                pvs_names_current_range,
            )[1]
        end
    end
    return res, output
end

function concatonate_i_true(fault_data, pvs_names_subset, range)
    i_true = []
    for (i, pvs_name) in enumerate(pvs_names_subset)
        i_true_fault = vcat(
            (fault_data[pvs_name][:ir_ground_truth])',
            (fault_data[pvs_name][:ii_ground_truth])',
        )
        i_true_fault = i_true_fault[:, range]
        if i == 1
            i_true = i_true_fault
        else
            i_true = hcat(i_true, i_true_fault)
        end
    end
    return i_true
end

function concatonate_pvs_names(pvs_names_subset, length_range)
    concatonated_pvs_names_list = []
    for (i, pvs_name) in enumerate(pvs_names_subset)
        if i == 1
            concatonated_pvs_names_list = fill(pvs_name, length_range)
        else
            concatonated_pvs_names_list =
                vcat(concatonated_pvs_names_list, fill(pvs_name, length_range))
        end
    end
    return concatonated_pvs_names_list
end

function concatonate_t(tsteps, pvs_names_subset, range)
    t = []
    for (i, pvs_name) in enumerate(pvs_names_subset)
        t_fault = (tsteps[range])'
        if i == 1
            t = t_fault
        else
            t = hcat(t, t_fault)
        end
    end
    return t
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
