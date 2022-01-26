using PowerSimulationNODE

params_data = NODETrainParams[]

push!(
    params_data,
    NODETrainParams(
        train_id = "1",
        maxiters = 2000,
        training_groups = [(tspan = (0.0, 1.0), multiple_shoot_group_size = 10, multiple_shoot_continuity_term = 100.0, batching_sample_factor = 1.0)],
        node_width = 5,
        node_layers = 2,
        node_activation = "relu",
        graphical_report_mode = 3,
    ),
)

push!(
    params_data,
    NODETrainParams(
        train_id = "2",
        maxiters = 2000,
        training_groups = [(tspan = (0.0, 1.0), multiple_shoot_group_size = 10, multiple_shoot_continuity_term = 100.0, batching_sample_factor = 1.0)],
        node_width = 10,
        node_layers = 2,
        node_activation = "relu",
        graphical_report_mode = 3,
    ),
)

push!(
    params_data,
    NODETrainParams(
        train_id = "3",
        maxiters = 2000,
        training_groups = [(tspan = (0.0, 1.0), multiple_shoot_group_size = 10, multiple_shoot_continuity_term = 100.0, batching_sample_factor = 1.0)],
        node_width = 15,
        node_layers = 2,
        node_activation = "relu",
        graphical_report_mode = 3,
    ),
)

push!(
    params_data,
    NODETrainParams(
        train_id = "4",
        maxiters = 2000,
        training_groups = [(tspan = (0.0, 1.0), multiple_shoot_group_size = 10, multiple_shoot_continuity_term = 100.0, batching_sample_factor = 1.0)],
        node_width = 20,
        node_layers = 2,
        node_activation = "relu",
        graphical_report_mode = 3,
    ),
)

push!(
    params_data,
    NODETrainParams(
        train_id = "5",
        maxiters = 2000,
        training_groups = [(tspan = (0.0, 1.0), multiple_shoot_group_size = 10, multiple_shoot_continuity_term = 100.0, batching_sample_factor = 1.0)],
        node_width = 5,
        node_layers = 3,
        node_activation = "relu",
        graphical_report_mode = 3,
    ),
)

push!(
    params_data,
    NODETrainParams(
        train_id = "6",
        maxiters = 2000,
        training_groups = [(tspan = (0.0, 1.0), multiple_shoot_group_size = 10, multiple_shoot_continuity_term = 100.0, batching_sample_factor = 1.0)],
        node_width = 10,
        node_layers = 3,
        node_activation = "relu",
        graphical_report_mode = 3,
    ),
)

push!(
    params_data,
    NODETrainParams(
        train_id = "7",
        maxiters = 2000,
        training_groups = [(tspan = (0.0, 1.0), multiple_shoot_group_size = 10, multiple_shoot_continuity_term = 100.0, batching_sample_factor = 1.0)],
        node_width = 15,
        node_layers = 3,
        node_activation = "relu",
        graphical_report_mode = 3,
    ),
)

push!(
    params_data,
    NODETrainParams(
        train_id = "8",
        maxiters = 2000,
        training_groups = [(tspan = (0.0, 1.0), multiple_shoot_group_size = 10, multiple_shoot_continuity_term = 100.0, batching_sample_factor = 1.0)],
        node_width = 20,
        node_layers = 3,
        node_activation = "relu",
        graphical_report_mode = 3,
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
    n_tasks = 6,
    force_generate_inputs = true,
)

generate_train_files(hpc_params)
run_parallel_train(hpc_params)
