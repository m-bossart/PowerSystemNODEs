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
    params_data = params_data,  #have the default time be infinite...
    project_folder = "PowerSystemNODEs",
    scratch_path = "/scratch/summit/mabo4366",
    n_tasks = length(params_data),
    time_limit = "23:59:59",        #add a time parameter which allocates (say 30 mins) to generating data, getting things instantiated, etc. 
    QoS = "normal",
    partition = "shas", #"shas-testing"
    force_generate_inputs = true,
    mb_per_cpu = 4800,
)

generate_train_files(hpc_params)    #when you generate the train data, first re-set the time limit in params before serializing it... 
##                                     #TODO - actually implement the time limit for training.
run_parallel_train(hpc_params)
