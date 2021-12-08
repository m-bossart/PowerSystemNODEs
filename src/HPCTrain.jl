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
#SBATCH --nodes={{n_nodes}}
#SBATCH --cpus-per-task=1
#
#SBATCH --time={{time_limit}}
#SBATCH --output={{{project_path}}}/job_output_%j.o
#SBATCH --error={{{project_path}}}/job_output_%j.e

# Check Dependencies
julia --project={{{project_path}}} -e 'using Pkg; Pkg.instantiate()'
julia --project={{{project_path}}} -e 'include("scripts/prepare_for_train.jl'

# Load Parallel
module load {{gnu_parallel_name}}

echo \$SLURM_JOB_NODELIST |sed s/\\,/\\\\n/g > hostfile

# --slf is needed to parallelize across all the cores on multiple nodes
parallel --jobs \$SLURM_CPUS_ON_NODE \\
    --slf hostfile \\
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
    n_nodes::Int
    params_data::Vector # TODO: return to Vector{NODETrainParams} after testing
    time_limit::String
    train_bash_file::String
end

function SavioHPCTrain(;
    username,
    params_data,
    project_folder = "PowerSystemNODEs",
    scratch_path = "/global/scratch/users",
    time_limit = "24:00:00",
    n_nodes = 1,
)
    return HPCTrain(
        username,
        "fc_emac",
        "savio_normal",
        "savio",
        project_folder,
        scratch_path,
        "gnu-parallel",
        n_nodes,
        params_data,
        time_limit,
        "",
    )
end

# Populated with info from: https://curc.readthedocs.io/en/latest/
function SummitHPCTrain(;
    username,
    params_data,
    project_folder = "PowerSystemNODEs",
    scratch_path = "/scratch/summit/",
    time_limit = "24:00:00",
    n_nodes = 1,
)
    return HPCTrain(
        username,
        "ucb-general", # The proper value is TBD
        "normal",
        "shas",
        project_folder,
        scratch_path,
        "gnu_parallel",
        n_nodes,
        params_data,
        time_limit,
        "",
    )
end

function generate_train_files(train::HPCTrain)
    data = Dict()
    data["username"] = train.username
    data["account"] = train.account
    data["QoS"] = train.QoS
    data["time_limit"] = train.time_limit
    data["partition"] = train.partition
    data["gnu_parallel_name"] = train.gnu_parallel_name
    data["project_path"] = joinpath(train.scratch_path, train.project_folder)
    data["n_nodes"] = train.n_nodes
    data["train_set_file"] =
        joinpath(train.scratch_path, train.project_folder, "train_files.lst")
    touch(data["train_set_file"])
    open(data["train_set_file"], "w") do file
        for param in train.params_data
            param_file_path = joinpath(
                train.scratch_path,
                train.project_folder,
                INPUT_FOLDER_NAME,
                "train_$(param.train_id).json",
            )
            touch(param_file_path)
            serialize(param, param_file_path)
            write(file, "$param_file_path\n")
        end
    end

    filename = HPC_TRAIN_FILE
    open(filename, "w") do io
        write(io, Mustache.render(bash_file_template, data))
    end
    return
end

function run_parallel_train(train::HPCTrain)
    return run(`sbatch $train.train_bash_file`)
end
