include("../src/SurrogateModels.jl")

optimizer_map = Dict("Adam" => ADAM, "Bfgs" => BFGS)

solver_map = Dict("Rodas4" => Rodas4)

sensealg_map = Dict("ForwardDiffSensitivity" => ForwardDiffSensitivity)

surr_map = Dict(
    "vsm_v_t_0" => vsm_v_t_0,
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

activation_map = Dict("relu" => relu)

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

    nn_input = 0
    nn_output = 2

    (inputs.node_inputs == "voltage") && (nn_input += 2)
    (inputs.node_feedback_current) && (nn_input += 2)
    nn_output += inputs.node_feedback_states
    nn_input += inputs.node_feedback_states
    Random.seed!(inputs.rng_seed)
    return build_nn(nn_input, nn_output, nn_width, nn_hidden, nn_activation)
end

function instantiate_M(inputs)
    n_differential = ODE_ORDER + 2 + inputs.node_feedback_states
    n_algebraic = 2

    return MassMatrix(n_differential, n_algebraic)
end

function instantiate_surr(surr, nn, Vm, Vθ)
    return (dx, x, p, t) -> surr(dx, x, p, t, nn, Vm, Vθ)
end

function instantiate_surr(inputs::NODETrainParams, nn, Vm, Vθ)
    if inputs.ode_model == "vsm"
        if inputs.node_inputs == "voltage"
            if inputs.node_feedback_current
                surr = surr_map[string("vsm_v_t_", inputs.node_feedback_states)]
                return instantiate_surr(surr, nn, Vm, Vθ)
            else
                return surr_map[string("vsm_v_f_", inputs.node_feedback_states)]
            end
        else
            @warn "node input type not found during surrogate instantiatiion"
        end
    else
        @warn "ode model not found during surrogate instantiatiion"
    end
end

#CLOSURE
function _loss_function(θ, y_actual, tsteps, weights, Ir_scale, Ii_scale, pred_function)
    y_predicted = pred_function(θ, tsteps)
    loss =
        (mae(y_predicted[1, :], y_actual[1, :]) / Ir_scale) +
        (mae(y_predicted[2, :], y_actual[2, :]) / Ii_scale) * weights[1] +
        (mse(y_predicted[1, :], y_actual[1, :]) / Ir_scale) +
        (mse(y_predicted[2, :], y_actual[2, :]) / Ii_scale) * weights[2]
    return loss, y_predicted
end

function instantiate_loss_function(weights, Ir_scale, Ii_scale, pred_function)
    return (θ, y_actual, tsteps) ->
        _loss_function(θ, y_actual, tsteps, weights, Ir_scale, Ii_scale, pred_function)
end

function _pred_function(θ, tsteps, p_fixed, solver, surr_prob, tols, sensealg, u₀)
    p = vcat(θ, p_fixed)
    _prob = remake(surr_prob, p = p, u0 = u₀)
    sol = solve(
        _prob,
        solver,
        abstol = tols[1],
        reltol = tols[2],
        saveat = tsteps,
        save_idxs = [I__IR_OUT, I__II_OUT, I__IR_FILTER, I__II_FILTER, I__IR_NN, I__II_NN], #first two for loss function, rest for data export
        sensealg = ForwardDiffSensitivity(),
    )
    return Array(sol)
end

function instantiate_pred_function(p_fixed, solver, surr_prob, tols, sensealg, u₀)
    return (θ, tsteps) ->
        _pred_function(θ, tsteps, p_fixed, solver, surr_prob, tols, sensealg, u₀)
end

function instantiate_cb!(output, lb_loss, exportmode, id, range_count)
    if exportmode == 3
        return (p, l, pred) -> _cb3!(p, l, pred, output, lb_loss, id, range_count)
    elseif exportmode == 2
        return (p, l, pred) -> _cb2!(p, l, pred, output, lb_loss, id, range_count)
    elseif exportmode == 1
        return (p, l, pred) -> _cb1!(p, l, pred, output, lb_loss, id, range_count)
    end
end

function _cb3!(p, l, pred, output, lb_loss, id, range_count)
    push!(output["loss"], (id, range_count, l))
    push!(output["parameters"], [p])
    push!(output["predictions"], (pred[1, :], pred[2, :]))
    output["total_iterations"] += 1
    @info "loss", l
    @info "p[end]", p[end]
    (l > lb_loss) && return false
    return true
end

function _cb2!(p, l, pred, output, lb_loss, id, range_count)
    push!(output["loss"], (id, range_count, l))
    output["total_iterations"] += 1
    @info "loss", l
    @info "p[end]", p[end]
    (l > lb_loss) && return false
    return true
end

function _cb1!(p, l, pred, output, lb_loss, id, range_count)
    output["total_iterations"] += 1
    @info "loss", l
    @info "p[end]", p[end]
    (l > lb_loss) && return false
    return true
end
