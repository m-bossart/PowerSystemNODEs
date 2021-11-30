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
#SBATCH --qos={{QoS}}(
#
# Partition:
#SBATCH --partition={{partition}}
#
#SBATCH --nodes={{n_nodes}}
#SBATCH --cpus-per-task=1
#
#SBATCH --output={{project_path}}/job_output_%j.o
#SBATCH --error={{project_path}}/job_output_%j.e

# Check Dependencies
julia --project={{project_path}} -e 'using Pkg; Pkg.instantiate("{{project_path}}")'

# Load Parallel
module load {{gnu_parallel_name}}

echo \$SLURM_JOB_NODELIST |sed s/\\,/\\\\n/g > hostfile

# --slf is needed to parallelize across all the cores on multiple nodes
parallel --jobs \$SLURM_CPUS_ON_NODE \\
    --slf hostfile \\
    --wd {{project_path}} \\
    --progress -a {{train_set_file}}.lst \\
    julia --project={{project_path}} {{project_path}}/scripts/train.jl {}
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
    params_data::Vector{NODETrainParams}
    input_data::NODETrainInputs
    train_bash_file::String
end

function SavioHPCTrain(;
    username,
    params_data,
    input_data = "train_data.json",
    project_folder = "PowerSystemNODEs",
    n_nodes = "1",
)
    return HPCTrain(
        username,
        "fc_emac",
        "savio_normal",
        "savio",
        project_folder,
        "/global/scratch/users",
        "gnu-parallel",
        n_nodes,
        params_data,
        input_data,
        "",
    )
end

# Populated with info from: https://curc.readthedocs.io/en/latest/
function SummitHPCTrain(;
    username,
    params_data,
    input_data = "train_data.json",
    project_folder = "PowerSystemNODEs",
    n_nodes = "1",
)
    return HPCTrain(
        username,
        "cu_allocation", # The proper value is TBD
        "normal",
        "shas",
        project_folder,
        "/scratch/summit/",
        "gnu_parallel",
        n_nodes,
        params_data,
        input_data,
        "",
    )
end

function generate_train_files(train::HPCTrain)
    data = Dict()
    data["username"] = train.username
    data["account"] = train.account
    data["QoS"] = train.QoS
    data["partition"] = train.partition
    data["gnu_parallel_name"] = train.gnu_parallel_name
    data["project_path"] = joinpath(train.scratch_path, train.project_folder)
    data["n_nodes"] = train.n_nodes
    data["train_set_file"] =
        joinpath(train.scratch_path, train.project_folder, "train_files.lst")
    open(data["train_set_file"], "w") do file
        for param in train.params_data
            param_file_path = joinpath(
                train.scratch_path,
                train.project_folder,
                INPUT_FOLDER_NAME,
                param.train_id,
                ".json",
            )
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
