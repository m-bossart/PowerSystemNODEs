using Pkg
Pkg.activate(".")
Pkg.instantiate()

using Mustache
include("../src/constants.jl")
include("../src/PowerSystemNODEs.jl")
include("../src/HPCTrain.jl")


params_data = NODETrainParams[]
for i in 1:2
    push!(params_data, NODETrainParams(train_id = string(i)))
end

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
)


generate_train_files(hpc_params)
run_parallel_train(hpc_params) 
