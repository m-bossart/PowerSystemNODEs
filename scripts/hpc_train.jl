using PowerSimulationNODE

params_data = NODETrainParams[]
no_change_params = Dict{Symbol, Any}()
change_params = Dict{Symbol, Any}()

#INDICATE CONSTANT, NON-DEFAULT PARAMETERS
no_change_params[:maxiters] = 10
no_change_params[:node_layers] = 1
no_change_params[:node_unobserved_states] = 8

#change_params[:sensealg] = ["ForwardDiff", "Zygote" ]
#change_params[:node_width] = [10, 200]
no_change_params[:sensealg] = "ForwardDiff"
no_change_params[:node_width] = 10

change_params[:solver] = ["Tsit5", "Rodas4", "TRBDF2"]
change_params[:solver_tols] = [(1e-5, 1e-2), (1e-7, 1e-4)]
change_params[:solver_sensealg] =
    ["InterpolatingAdjoint", "InterpolatingAdjoint_checkpointing"]

#INDICATE PARAMETES TO ITERATE OVER COMBINATORIALLY 

#SPECIAL HANDLING TO BUILD ITERATOR FOR TRAINING GROUPS 
no_change_fields = Dict{Symbol, Any}()
change_fields = Dict{Symbol, Any}()
no_change_fields[:tspan] = (0.0, 1.0)
no_change_fields[:training_groups] = 1
no_change_fields[:multiple_shoot_continuity_term] = (0.0, 1.0)
no_change_fields[:batching_sample_factor] = 0.5
change_fields[:shoot_times] = [[], [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]]
change_params[:training_groups] =
    build_training_groups_list(no_change_fields, change_fields)

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
    n_tasks = length(params_data),
    time_limit = "24:00:00",
    QoS = "normal",
    partition = "shas", #"shas-testing"
    force_generate_inputs = true,
    mb_per_cpu = 4800,
)

generate_train_files(hpc_params)
##
run_parallel_train(hpc_params)
