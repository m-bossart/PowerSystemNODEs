#using Pkg
#Pkg.activate(".")

using Mustache
include("../src/constants.jl")
include("../src/HPCTrain.jl")

struct TestParams
    train_id::String
end

test = [TestParams("test"), TestParams("test2")]

hpc_params = SavioHPCTrain(;
    username = "jdlara",
    params_data = test,
    project_folder = "test",
    scratch_path = "/global/home/users/jdlara",
)

generate_train_files(hpc_params)
