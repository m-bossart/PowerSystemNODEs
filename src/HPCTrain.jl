using Mustache

const bash_file_template = """
#!/bin/bash
# Job name:
#SBATCH --job-name=NODE_train
#
# Account:
#SBATCH --account={{account}}
#
# QoS:
#SBATCH --qos={{QoS}}
#
# Partition:
#SBATCH --partition={{partition}}
#
# Number of MPI tasks requested:
#SBATCH --ntasks={{n_tasks}} \n{{#n_nodes}}#\n#SBATCH --nodes={{n_nodes}} \n {{/n_nodes}}
#
# Processors per task (for future parallel training code):
#SBATCH --cpus-per-task={{n_cpus_per_task}}
#
#SBATCH --time={{time_limit}}
#SBATCH --output={{{project_path}}}/job_output_%j.o
#SBATCH --error={{{project_path}}}/job_output_%j.e

# Check Dependencies
julia --project={{{project_path}}} -e 'using Pkg; Pkg.instantiate()'
{{#force_generate_inputs}} julia --project={{{project_path}}} {{{project_path}}}/scripts/prepare_for_train.jl "true" {{/force_generate_inputs}}

# Load Parallel
module load {{gnu_parallel_name}}

{{#n_nodes}}
# --slf is needed to parallelize across all the cores on multiple nodes
dataecho \$SLURM_JOB_NODELIST |sed s/\\,/\\\\n/g > hostfile
{{/n_nodes}}

{{#n_nodes}}parallel --jobs \$SLURM_CPUS_ON_NODE --slf hostfile \\ {{/n_nodes}}{{^n_nodes}}parallel --jobs \$SLURM_NPROCS \\{{/n_nodes}}
    --wd {{{project_path}}} \\
    --progress -a {{{train_set_file}}}\\
    --joblog {{{project_path}}}/hpc_train.log \\
    julia --project={{{project_path}}} {{{project_path}}}/scripts/train_node.jl {}
"""

struct HPCTrain
    username::String
    account::String
    QoS::String
    partition::String
    project_folder::String
    # TODO: Coordinate properly with the data in the inputs vector base_path field
    scratch_path::String
    gnu_parallel_name::String
    n_tasks::Int
    n_nodes::Union{Int, Nothing}
    n_cpus_per_task::Int
    params_data::Vector # TODO: return to Vector{NODETrainParams} after testing
    time_limit::String
    train_bash_file::String
    force_generate_inputs::Bool
end

function SavioHPCTrain(;
    username,
    params_data,
    project_folder = "PowerSystemNODEs",
    scratch_path = "/global/scratch/users",
    time_limit = "24:00:00",
    n_tasks = 1,
    force_generate_inputs = false,
    n_nodes = nothing, # Use with caution in Savio, it can lead to over subscription of nodes
)

    # Default until we parallelize training code
    n_cpus_per_task = 1
    return HPCTrain(
        username,
        "fc_emac",
        "savio_normal",
        "savio",
        project_folder,
        scratch_path,
        "gnu-parallel",
        n_tasks,
        n_nodes,
        n_cpus_per_task,
        params_data,
        time_limit,
        joinpath(scratch_path, project_folder, HPC_TRAIN_FILE),
        force_generate_inputs,
    )
end

# Populated with info from: https://curc.readthedocs.io/en/latest/
function SummitHPCTrain(;
    username,
    params_data,
    project_folder = "PowerSystemNODEs",
    scratch_path = "/scratch/summit/",
    time_limit = "24:00:00",
    n_tasks = 1,  #default to parallelize across all tasks 
    force_generate_inputs = false,
)
    # Default until we parallelize training code
    n_cpus_per_task = 1
    return HPCTrain(
        username,
        "ucb-general", # Get allocation
        "normal",
        "shas",
        project_folder,
        scratch_path,
        "gnu_parallel",
        n_tasks,
        nothing, # Default to nothing on Summit since it doesn't dispatch on ssh login
        n_cpus_per_task,
        params_data,
        time_limit,
        joinpath(scratch_path, project_folder, HPC_TRAIN_FILE),
        force_generate_inputs,
    )
end

function generate_train_files(train::HPCTrain)
    mkpath(INPUT_FOLDER_NAME)
    mkpath(INPUT_SYSTEM_FOLDER_NAME)
    mkpath(OUTPUT_FOLDER_NAME)
    touch(HPC_TRAIN_FILE)
    #touch(joinpath(INPUT_FOLDER_NAME,"data.json"))

    data = Dict()
    data["username"] = train.username
    data["account"] = train.account
    data["QoS"] = train.QoS
    data["time_limit"] = train.time_limit
    data["partition"] = train.partition
    data["gnu_parallel_name"] = train.gnu_parallel_name
    data["project_path"] = joinpath(train.scratch_path, train.project_folder)
    data["n_tasks"] = train.n_tasks
    data["n_cpus_per_task"] = train.n_cpus_per_task

    data["n_nodes"] = train.n_nodes
    if !isnothing(train.n_nodes)
        data["n_nodes"] = train.n_nodes
    end

    data["force_generate_inputs"] = train.force_generate_inputs
    train_set_folder = joinpath(train.scratch_path, train.project_folder)
    if !ispath(train_set_folder)
        mkpath(train_set_folder)
    end
    data["train_set_file"] = joinpath(train_set_folder, "train_files.lst")

    open(data["train_set_file"], "w") do file
        for param in train.params_data
            param_file_path = joinpath(
                train.scratch_path,
                train.project_folder,
                INPUT_FOLDER_NAME,
                "train_$(param.train_id).json",
            )
            if !isfile(param_file_path)
                touch(param_file_path)
            end
            serialize(param, param_file_path)
            write(file, "$param_file_path\n")
        end
    end

    data["force_generate_inputs"] = train.force_generate_inputs
    open(train.train_bash_file, "w") do io
        write(io, Mustache.render(bash_file_template, data))
    end
    return
end

function run_parallel_train(train::HPCTrain)
    bash_file = train.train_bash_file
    return run(`sbatch $bash_file`)
end
