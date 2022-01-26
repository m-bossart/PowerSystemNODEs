using PowerSimulationNODE

params_data = NODETrainParams[]

push!(
    params_data,
    NODETrainParams()
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
    n_tasks = 1,
    force_generate_inputs = true,
)

generate_train_files(hpc_params)
run_parallel_train(hpc_params)
