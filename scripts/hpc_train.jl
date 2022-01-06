using Mustache

include("../src/PowerSystemNODEs.jl")
include("../src/constants.jl")
include("../src/HPCTrain.jl")

params_data = NODETrainParams[]

push!(
    params_data,
    NODETrainParams(
        ode_model = "none",
        train_id = "1",
        node_width = 2,
        node_layers = 2,
        node_feedback_states = 0,
    ),
)
push!(
    params_data,
    NODETrainParams(
        ode_model = "none",
        train_id = "2",
        node_width = 3,
        node_layers = 2,
        node_feedback_states = 0,
    ),
)
push!(
    params_data,
    NODETrainParams(
        ode_model = "none",
        train_id = "3",
        node_width = 4,
        node_layers = 2,
        node_feedback_states = 0,
    ),
)
push!(
    params_data,
    NODETrainParams(
        ode_model = "none",
        train_id = "4",
        node_width = 5,
        node_layers = 2,
        node_feedback_states = 0,
    ),
)
push!(
    params_data,
    NODETrainParams(
        ode_model = "none",
        train_id = "5",
        node_width = 2,
        node_layers = 3,
        node_feedback_states = 0,
    ),
)
push!(
    params_data,
    NODETrainParams(
        ode_model = "none",
        train_id = "6",
        node_width = 2,
        node_layers = 4,
        node_feedback_states = 0,
    ),
)
push!(
    params_data,
    NODETrainParams(
        ode_model = "none",
        train_id = "7",
        node_width = 2,
        node_layers = 5,
        node_feedback_states = 0,
    ),
)
push!(
    params_data,
    NODETrainParams(
        ode_model = "none",
        train_id = "8",
        node_width = 2,
        node_layers = 2,
        node_feedback_states = 1,
    ),
)
push!(
    params_data,
    NODETrainParams(
        ode_model = "none",
        train_id = "9",
        node_width = 2,
        node_layers = 2,
        node_feedback_states = 2,
    ),
)
push!(
    params_data,
    NODETrainParams(
        ode_model = "none",
        train_id = "10",
        node_width = 2,
        node_layers = 2,
        node_feedback_states = 3,
    ),
)

#=
 hpc_params = SavioHPCTrain(;
    username = "jdlara",
    params_data = params_data,
    project_folder = "PowerSystemNODEs",
    scratch_path = "/global/home/users/jdlara",
)
  =#

hpc_params = SummitHPCTrain(;
    username = "mabo4366",
    params_data = params_data,
    project_folder = "PowerSystemNODEs",
    scratch_path = "/scratch/summit/mabo4366",
    n_tasks = 10,
    force_generate_inputs = true, 
)

generate_train_files(hpc_params)
run_parallel_train(hpc_params)
