const optimizer_map = Dict("Adam" => ADAM, "Bfgs" => BFGS)

const solver_map = Dict("Rodas4" => Rodas4)

const sensealg_map = Dict("ForwardDiffSensitivity" => ForwardDiffSensitivity)

const surr_map = Dict(
    "vsm_v_t_0" => vsm_v_t_0,
    "none_v_t_0" => none_v_t_0,
    "none_v_t_1" => none_v_t_1,
    "none_v_t_2" => none_v_t_2,
    "none_v_t_3" => none_v_t_3,
    "none_v_t_4" => none_v_t_4,
    "none_v_t_5" => none_v_t_5,

    #=  "vsm_v_t_3" => vsm_v_t_3,
        "vsm_v_t_4" => vsm_v_t_4,
        "vsm_v_t_5" => vsm_v_t_5,
        "vsm_v_t_6" => vsm_v_t_6,
        "vsm_v_t_7" => vsm_v_t_7,
        "vsm_v_t_8" => vsm_v_t_8,
        "vsm_v_t_9" => vsm_v_t_9,
        "vsm_v_t_10" => vsm_v_t_10,
        "vsm_v_t_11" => vsm_v_t_11,
        "vsm_v_t_12" => vsm_v_t_12,
        "vsm_v_t_13" => vsm_v_t_13,
        "vsm_v_t_14" => vsm_v_t_14,
        "vsm_v_t_15" => vsm_v_t_15,
        "vsm_v_t_16" => vsm_v_t_16,
        "vsm_v_t_17" => vsm_v_t_17,
        "vsm_v_t_18" => vsm_v_t_18,
        "vsm_v_t_19" => vsm_v_t_19, =#
)

const activation_map = Dict("relu" => relu)

function instantiate_solver(inputs)
    return solver_map[inputs.solver]()
end

function instantiate_sensealg(inputs)
    return sensealg_map[inputs.sensealg]()
end

function instantiate_optimizer(inputs)
    if inputs.optimizer == "Adam"
        return optimizer_map[inputs.optimizer](inputs.optimizer_η)
    elseif inputs.optimizer == "Bfgs"
        return optimizer_map[inputs.optimizer]()
    end
end

function instantiate_optimizer_adjust(inputs)
    if inputs.optimizer_adjust == "Adam"
        return optimizer_map[inputs.optimizer_adjust](inputs.optimizer_adjust_η)
    elseif inputs.optimizer_adjust == "Bfgs"
        return optimizer_map[inputs.optimizer_adjust]()
    end
end

function instantiate_nn(inputs)
    nn_activation = activation_map[inputs.node_activation]
    nn_hidden = inputs.node_layers
    nn_width = inputs.node_width

    nn_input = 4   #P_pf, Q_pf, V_pf, θ_pf
    nn_output = 2

    (inputs.node_inputs == "voltage") && (nn_input += 2)
    (inputs.node_feedback_current) && (nn_input += 2)
    nn_output += inputs.node_feedback_states
    nn_input += inputs.node_feedback_states
    Random.seed!(inputs.rng_seed)
    @warn "NN size parameters" nn_input, nn_output, nn_width, nn_hidden
    return build_nn(nn_input, nn_output, nn_width, nn_hidden, nn_activation)
end

function instantiate_M(inputs)
    if inputs.ode_model == "vsm"
        ODE_ORDER = 19
        N_ALGEBRAIC = 2
    elseif inputs.ode_model == "none"
        ODE_ORDER = 0
        N_ALGEBRAIC = 0
    else
        @error "ODE order unknown for ODE model provided"
    end
    n_differential = ODE_ORDER + 2 + inputs.node_feedback_states

    return MassMatrix(n_differential, N_ALGEBRAIC)
end

function instantiate_surr(surr, nn, Vm, Vθ)
    return (dx, x, p, t) -> surr(dx, x, p, t, nn, Vm, Vθ)
end

function instantiate_surr(inputs::NODETrainParams, nn, Vm, Vθ)
    if inputs.ode_model == "vsm"
        N_ALGEBRAIC_STATES = 2
        ODE_ORDER = 19
        if inputs.node_inputs == "voltage"
            if inputs.node_feedback_current
                surr = surr_map[string("vsm_v_t_", inputs.node_feedback_states)]
                return instantiate_surr(surr, nn, Vm, Vθ), N_ALGEBRAIC_STATES, ODE_ORDER
            else
                surr = surr_map[string("vsm_v_f_", inputs.node_feedback_states)]
                return instantiate_surr(surr, nn, Vm, Vθ), N_ALGEBRAIC_STATES, ODE_ORDER
            end
        else
            @warn "node input type not found during surrogate instantiatiion"
        end
    elseif inputs.ode_model == "none"
        N_ALGEBRAIC_STATES = 0
        ODE_ORDER = 0
        if inputs.node_inputs == "voltage"
            if inputs.node_feedback_current
                surr = surr_map[string("none_v_t_", inputs.node_feedback_states)]
                return instantiate_surr(surr, nn, Vm, Vθ), N_ALGEBRAIC_STATES, ODE_ORDER
            else
                surr = surr_map[string("none_v_f_", inputs.node_feedback_states)]
                return instantiate_surr(surr, nn, Vm, Vθ), N_ALGEBRAIC_STATES, ODE_ORDER
            end
        else
            @warn "node input type not found during surrogate instantiatiion"
        end
    else
        @warn "ode model not found during surrogate instantiatiion"
    end
end

function _loss_function(
    θ,
    y_actual,
    tsteps,
    weights,
    Ir_scale,
    Ii_scale,
    pred_function,
    pvs_names,
)
    y_predicted = pred_function(θ, tsteps, pvs_names)     #Careful of dict ordering?                  

    if size(y_predicted) == size(y_actual)
        loss =
            (mae(y_predicted[1, :], y_actual[1, :]) / Ir_scale) +
            (mae(y_predicted[2, :], y_actual[2, :]) / Ii_scale) * weights[1] +
            (mse(y_predicted[1, :], y_actual[1, :]) / Ir_scale) +
            (mse(y_predicted[2, :], y_actual[2, :]) / Ii_scale) * weights[2]
    else
        loss = Inf
        @warn "Unstable run detected, assigning infinite loss"
    end

    return loss, y_predicted
end

function instantiate_loss_function(weights, Ir_scale, Ii_scale, pred_function)
    return (θ, y_actual, tsteps, pvs_names) -> _loss_function(
        θ,
        y_actual,
        tsteps,
        weights,
        Ir_scale,
        Ii_scale,
        pred_function,
        pvs_names,
    )
end

function _pred_function(θ, tsteps, P, solver, surr_prob, tols, sensealg, u₀)
    P.nn = θ
    p = vectorize(P)
    _prob = remake(surr_prob, p = p, u0 = eltype(p).(u₀))
    sol = solve(    #BUG - dimension mismatch 2-4. 
        _prob,
        solver,
        abstol = tols[1],
        reltol = tols[2],
        saveat = tsteps,
        save_idxs = [1, 2],  #First 2 states always the output currents?  , I__II_OUT, I__IR_FILTER, I__II_FILTER, I__IR_NN, I__II_NN], #first two for loss function, rest for data export TODO - should not be constant, depends on surrogate model 
        sensealg = ForwardDiffSensitivity(),
    )
    #@warn "sol", Array(sol)
    return Array(sol)
end

function full_array_pred_function(
    θ,
    tsteps,
    solver,
    pvs_names_subset,
    fault_data,
    tols,
    sensealg,
    pvs_names,
)
    full_array = []
    for (i, pvs_name) in enumerate(pvs_names_subset)
        surr_prob = fault_data[pvs_name][:surr_problem]
        u₀ = surr_prob.u0
        P = fault_data[pvs_name][:P]
        selector = [pvs_name .== name for name in pvs_names]
        t_steps_subset = tsteps[selector]
        if i == 1
            full_array =
                _pred_function(θ, t_steps_subset, P, solver, surr_prob, tols, sensealg, u₀)
        else
            full_array = hcat(
                full_array,
                _pred_function(θ, t_steps_subset, P, solver, surr_prob, tols, sensealg, u₀),
            )
        end
    end
    return full_array
end

function instantiate_pred_function(solver, pvs_names_subset, fault_data, tols, sensealg)
    return (θ, tsteps, pvs_names) -> full_array_pred_function(
        θ,
        tsteps,
        solver,
        pvs_names_subset,
        fault_data,
        tols,
        sensealg,
        pvs_names,
    )
end

function instantiate_cb!(output, lb_loss, exportmode, range_count, pvs_names, t_prediction)
    if exportmode == 3
        return (p, l, pred) ->
            _cb3!(p, l, pred, output, lb_loss, range_count, pvs_names, t_prediction)
    elseif exportmode == 2
        return (p, l, pred) ->
            _cb2!(p, l, pred, output, lb_loss, range_count, pvs_names, t_prediction)
    elseif exportmode == 1
        return (p, l, pred) ->
            _cb1!(p, l, pred, output, lb_loss, range_count, pvs_names, t_prediction)
    end
end

function _cb3!(p, l, pred, output, lb_loss, range_count, pvs_names, t_prediction)
    push!(output["loss"], (collect(pvs_names), range_count, l))
    push!(output["parameters"], [p])
    push!(output["predictions"], (vec(Float64.(t_prediction)), pred[1, :], pred[2, :]))
    output["total_iterations"] += 1
    @info "loss", l
    @info "p[end]", p[end]
    (l > lb_loss) && return false
    return true
end

function _cb2!(p, l, pred, output, lb_loss, range_count, pvs_names, t_prediction)
    push!(output["loss"], (collect(pvs_names), range_count, l))
    output["total_iterations"] += 1
    @info "loss", l
    @info "p[end]", p[end]
    (l > lb_loss) && return false
    return true
end

function _cb1!(p, l, pred, output, lb_loss, range_count, pvs_names, t_prediction)
    output["total_iterations"] += 1
    @info "loss", l
    @info "p[end]", p[end]
    (l > lb_loss) && return false
    return true
end
