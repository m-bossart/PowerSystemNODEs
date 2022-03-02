using PowerSimulationNODE
using Base.Cartesian
##
function build_params_list!(params_data, no_change_params, change_params)
    train_id = 1
    starting_dict = no_change_params
    dims = []
    for (k, v) in change_params
        push!(dims, length(v))
    end
    dims_tuple = tuple(dims...)
    iterator = CartesianIndices(dims_tuple)
    for i in iterator
        for (j, (key, value)) in enumerate(change_params)
            starting_dict[key] = value[i[j]]
        end
        starting_dict[:train_id] = string(train_id)
        push!(params_data, NODETrainParams(; starting_dict...))
        starting_dict = no_change_params
        train_id += 1
    end
end

function build_training_group(training_group_dict)
    training_group = []
    for i in 1:training_group_dict[:training_groups]
        tspan = (training_group_dict[:tspan][1], training_group_dict[:tspan][2] / i)
        shoot_times = filter(x -> x < tspan[2], training_group_dict[:shoot_times])
        multiple_shoot_continuity_term =
            training_group_dict[:multiple_shoot_continuity_term]
        batching_sample_factor = training_group_dict[:batching_sample_factor]
        push!(
            training_group,
            (
                tspan = tspan,
                shoot_times = shoot_times,
                multiple_shoot_continuity_term = multiple_shoot_continuity_term,
                batching_sample_factor = batching_sample_factor,
            ),
        )
    end
    return training_group
end
function build_training_groups_list(no_change_fields, change_fields)
    training_groups_list = []
    starting_dict = no_change_fields
    dims = []
    for (k, v) in change_fields
        push!(dims, length(v))
    end
    dims_tuple = tuple(dims...)
    iterator = CartesianIndices(dims_tuple)
    for i in iterator
        for (j, (key, value)) in enumerate(change_fields)
            starting_dict[key] = value[i[j]]
        end
        push!(training_groups_list, build_training_group(starting_dict))
    end

    return training_groups_list
end

params_data = NODETrainParams[]
no_change_params = Dict{Symbol, Any}()
change_params = Dict{Symbol, Any}()

#INDICATE CONSTANT, NON-DEFAULT PARAMETERS
no_change_params[:maxiters] = 2000
no_change_params[:node_layers] = 3
no_change_params[:node_width] = 8

#INDICATE PARAMETES TO ITERATE OVER COMBINATORIALLY 
change_params[:node_unobserved_states] = [0, 8]
change_params[:optimizer_η] = [0.001, 0.01]
change_params[:learn_initial_condition_unobserved_states] = [true, false]

#SPECIAL HANDLING TO BUILD ITERATOR FOR TRAINING GROUPS 
no_change_fields = Dict{Symbol, Any}()
change_fields = Dict{Symbol, Any}()
no_change_fields[:tspan] = (0.0, 1.0)
change_fields[:training_groups] = [1, 2]
change_fields[:shoot_times] = [
    [],
    [0.2, 0.4, 0.6, 0.8],
    [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
    [
        0.05,
        0.1,
        0.15,
        0.2,
        0.25,
        0.3,
        0.35,
        0.4,
        0.45,
        0.5,
        0.55,
        0.6,
        0.65,
        0.7,
        0.75,
        0.8,
        0.85,
        0.9,
        0.95,
    ],
]
change_fields[:multiple_shoot_continuity_term] = [10, 100]
change_fields[:batching_sample_factor] = [0.5, 1.0]
change_params[:training_groups] =
    build_training_groups_list(no_change_fields, change_fields)

build_params_list!(params_data, no_change_params, change_params)
@warn length(params_data)
##
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
    scratch_path = pwd(), #"/scratch/summit/mabo4366",
    n_tasks = length(params_data),
    time_limit = "24:00:00",
    QoS = "normal",
    partition = "shas", #"shas-testing"
    force_generate_inputs = true,
    mb_per_cpu = 4800,
)
##
generate_train_files(hpc_params)
run_parallel_train(hpc_params)
