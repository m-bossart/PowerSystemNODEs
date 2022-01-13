using PowerSimulationNODE

params_data = NODETrainParams[]

push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 10,
        ode_model = "none",
        train_id = "1",
        node_width = 2,
        node_layers = 2,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 10,
        ode_model = "none",
        train_id = "2",
        node_width = 3,
        node_layers = 2,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 10,
        ode_model = "none",
        train_id = "3",
        node_width = 4,
        node_layers = 2,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 10,
        ode_model = "none",
        train_id = "4",
        node_width = 5,
        node_layers = 2,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 10,
        ode_model = "none",
        train_id = "5",
        node_width = 2,
        node_layers = 3,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 10,
        ode_model = "none",
        train_id = "6",
        node_width = 2,
        node_layers = 4,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 10,
        ode_model = "none",
        train_id = "7",
        node_width = 2,
        node_layers = 5,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 10,
        ode_model = "none",
        train_id = "8",
        node_width = 2,
        node_layers = 2,
        node_feedback_states = 1,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 10,
        ode_model = "none",
        train_id = "9",
        node_width = 2,
        node_layers = 2,
        node_feedback_states = 2,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 10,
        ode_model = "none",
        train_id = "10",
        node_width = 2,
        node_layers = 2,
        node_feedback_states = 3,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 10,
        ode_model = "none",
        train_id = "11",
        node_width = 5,
        node_layers = 5,
        node_feedback_states = 3,
        graphical_report_mode = 3,
    ),
)

push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 5,
        ode_model = "none",
        train_id = "12",
        node_width = 2,
        node_layers = 2,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 5,
        ode_model = "none",
        train_id = "13",
        node_width = 3,
        node_layers = 2,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 5,
        ode_model = "none",
        train_id = "14",
        node_width = 4,
        node_layers = 2,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 5,
        ode_model = "none",
        train_id = "15",
        node_width = 5,
        node_layers = 2,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 5,
        ode_model = "none",
        train_id = "16",
        node_width = 2,
        node_layers = 3,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 5,
        ode_model = "none",
        train_id = "17",
        node_width = 2,
        node_layers = 4,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 5,
        ode_model = "none",
        train_id = "18",
        node_width = 2,
        node_layers = 5,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 5,
        ode_model = "none",
        train_id = "19",
        node_width = 2,
        node_layers = 2,
        node_feedback_states = 1,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 5,
        ode_model = "none",
        train_id = "20",
        node_width = 2,
        node_layers = 2,
        node_feedback_states = 2,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 5,
        ode_model = "none",
        train_id = "21",
        node_width = 2,
        node_layers = 2,
        node_feedback_states = 3,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 5,
        ode_model = "none",
        train_id = "22",
        node_width = 5,
        node_layers = 5,
        node_feedback_states = 3,
        graphical_report_mode = 3,
    ),
)

push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 2 ,
        ode_model = "none",
        train_id = "23",
        node_width = 2,
        node_layers = 2,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 2 ,
        ode_model = "none",
        train_id = "24",
        node_width = 3,
        node_layers = 2,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 2 ,
        ode_model = "none",
        train_id = "25",
        node_width = 4,
        node_layers = 2,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 2 ,
        ode_model = "none",
        train_id = "26",
        node_width = 5,
        node_layers = 2,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 2 ,
        ode_model = "none",
        train_id = "27",
        node_width = 2,
        node_layers = 3,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 2 ,
        ode_model = "none",
        train_id = "28",
        node_width = 2,
        node_layers = 4,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 2 ,
        ode_model = "none",
        train_id = "29",
        node_width = 2,
        node_layers = 5,
        node_feedback_states = 0,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 2 ,
        ode_model = "none",
        train_id = "30",
        node_width = 2,
        node_layers = 2,
        node_feedback_states = 1,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 2 ,
        ode_model = "none",
        train_id = "31",
        node_width = 2,
        node_layers = 2,
        node_feedback_states = 2,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 2 ,
        ode_model = "none",
        train_id = "32",
        node_width = 2,
        node_layers = 2,
        node_feedback_states = 3,
        graphical_report_mode = 3,
    ),
)
push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 2 ,
        ode_model = "none",
        train_id = "33",
        node_width = 5,
        node_layers = 5,
        node_feedback_states = 3,
        graphical_report_mode = 3,
    ),
)

push!(
    params_data,
    NODETrainParams(
        maxiters = 20000,
        groupsize_steps = 2 ,
        ode_model = "none",
        train_id = "33_repeat",
        node_width = 5,
        node_layers = 5,
        node_feedback_states = 3,
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
    n_tasks = 10,
    force_generate_inputs = true,
)

generate_train_files(hpc_params)
run_parallel_train(hpc_params)
