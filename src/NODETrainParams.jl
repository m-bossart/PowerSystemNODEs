"""
# Fields
- `train_id::Int64`: id for the training instance, used for naming output data folder.
- `solver::String`: solver used for the NODE problem. Valid Values ["Rodas4"]
- `solver_tols:: Tuple{Float64, Float64}`: solver tolerances (abstol, reltol).
- `sensealg::String`: sensitivity algorithm used in training. Valid Values ["ForwardDiffSensitivity"]
- `optimizer::["Adam", "Bfgs"]`: main optimizer used in training.
- `optimizer_η::Float64`: Learning rate for Adam (amount by which gradients are discounted before updating weights). Ignored if Adam is not the optimizer.
- `optimizer_adjust::String: optimizer used for final adjustments (2nd stage). Valid values ["Adam", "Bfgs", "nothing"].
- `optimizer_adjust_η::Float64`: Learning rate for Adam (amount by which gradients are discounted before updating weights). Ignored if Adam is not the optimizer.
- `maxiters::Int64`: The maximum possible iterations for the entire training instance. If `lb_loss = 0` and `optimizer = "Adam"` the training should never exit early and maxiters will be hit.
    Note that the number of saved data points can exceed maxiters because there is an additional callback at the end of each individual optimization.
- `lb_loss::Float64`: If the value of the loss function moves below lb_loss during training, the current optimization ends (current range).
- `batching::Bool`: If `batching = false` the full set of data points are used for each training step.
- `batching_factor::Float64`: The number of data points in the current range is multiplied by `batching_factor` to get the size of the batch. Batches of this size are used sequentially in time.
    The final batch is used even if it is incomplete.
- `rng_seed::Int64`: Seed for the random number generator used for initializing the NN for reproducibility across training runs.
- `groupsize_steps::Int64`: Number of data-points in each extension of the range of data used.
- `groupsize_faults::Int64`: Number of faults trained on simultaneous `1`:sequential training. if equal to number of pvs in sys_train, parallel training.
- `loss_function_weights::Tuple{Float64, Float64}`: weights used for loss function `(mae_weight, mse_weight)`.
- `loss_function_scale::String`: Scaling of the loss function.  `"range"`: the range of the real current and imaginary current are used to scale both the mae. Valid values ["range", "none"]
    and mse portions of the loss function. The goal is to give equal weight to real and imaginary components even if the magnitude of the disturbance differs. `"none"`: no additional scaling applied.
- `ode_model::String ["none","vsm"]`: The ode model used in conjunction with the NODE during training. `"none"` uses a purely data driven NODE surrogate. 
- `node_input_scale::Float64`: Scale factor on the voltage input to the NODE. Does not apply to other inputs (ie the feedback states).
- `node_output_scale::Float64`: Scale factor on the current output of the NODE. Does not apply to other outputs (ie the feedback states).
- `node_inputs::["voltage"]`: Determines the physical states which are inputs to the NODE. Ideally, only voltage to remain as general as possible.
- `node_feedback_states::Int64`: Number of feedback states in the NODE. Does not include the output current states which can be feedback if `node_feedback_current = true`.
- `node_feedback_current::Bool`: Determines if current is also a feedback state.
- `node_layers::Int64`: Number of hidden layers in the NODE. Does not include the input or output layer.
- `node_width::Int64`: Number of neurons in each hidden layer. Each hidden layer has the same number of neurons. The width of the input and output layers are determined by the combination of other parameters.
- `node_activation::String`: Activation function for NODE. The output layer always uses the identity activation. Valid Values ["relu"]
- `output_mode::Int`: `1`: do not collect any data during training, only save high-level data related to training and final results `2`: Same as `1`, also save value of loss throughout training. Valid values [1,2,3]
    `3`: same as `2`, also save parameters and predictions during training.
- `base_path:String`: Directory for training where input data is found and output data is written.
- `input_data_path:String`: From `base_path`, the directory for input data.
- `output_data_path:String`: From `base_path`, the directory for saving output data.
- `verify_psid_node_off:Bool`: `true`: before training, check that the surrogate with NODE turned off matches the data provided from PSID simulation.
- `graphical_report_mode:Int64`: `0`: do not generate plots. `1`: plot final result only. `2` plot for transitions between faults. `3`: plot for transitions between ranges. `4`: plot for every train iteration.
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
    graphical_report_mode::Int64
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
    graphical_report_mode = 0,
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
        graphical_report_mode,
    )
end

function NODETrainParams(file::AbstractString)
    return JSON3.read(read(file), NODETrainParams)
end

function read_input_data(pvs, d)
    id = get_name(pvs)
    tsteps = Float64.(d[:tsteps])
    i_ground_truth = vcat(Float64.(d[:ir_ground_truth])', Float64.(d[:ii_ground_truth])')
    i_node_off = vcat(Float64.(d[:ir_node_off])', Float64.(d[:ii_node_off])')
    p_ode = Float64.(d[:p_ode])
    x₀ = Float64.(d[:x₀])
    p_V₀ = Float64.(d[:V₀])
    return id, tsteps, i_ground_truth, i_node_off, p_ode, x₀, p_V₀
end

function serialize(inputs::NODETrainParams, file_path::String)
    open(file_path, "w") do io
        JSON3.write(io, inputs)
    end
    return
end
