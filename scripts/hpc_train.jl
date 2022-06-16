using PowerSimulationNODE

params_data = TrainParams[]
no_change_params = Dict{Symbol, Any}()
change_params = Dict{Symbol, Any}()

#INDICATE CONSTANT, NON-DEFAULT PARAMETERS
no_change_params[:maxiters] = 100

#INDICATE PARAMETES TO ITERATE OVER COMBINATORIALLY 
change_params[:optimizer] = [
    (
        sensealg = "Zygote",
        primary = "Adam",
        primary_η = 0.0001,
        adjust = "nothing",
        adjust_η = 0.0,
    ),
    (
        sensealg = "Zygote",
        primary = "Adam",
        primary_η = 0.0001,
        adjust = "nothing",
        adjust_η = 0.0,
    ),
    (
        sensealg = "Zygote",
        primary = "Adam",
        primary_η = 0.0001,
        adjust = "nothing",
        adjust_η = 0.0,
    ),
]

build_params_list!(params_data, no_change_params, change_params)
@warn "Number of trainings:", length(params_data)

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
    time_limit_train = "23:59:59",       
    time_limit_generate_data = "00:30:00", 
    QoS = "normal",
    partition = "shas", #"shas-testing"
    force_generate_inputs = true,
    mb_per_cpu = 4800,
)

generate_train_files(hpc_params)   
##                                   
run_parallel_train(hpc_params)
