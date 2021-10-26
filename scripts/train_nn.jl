using Pkg
Pkg.activate(".")
using Revise
using DifferentialEquations
using DiffEqSensitivity
using Logging
using PowerSystems
using PowerSimulationsDynamics
using GalacticOptim
using Plots
using IterTools
using NLsolve
using DiffEqFlux: group_ranges
using DiffEqFlux
using Flux
using Flux.Losses: mae, mse
using ForwardDiff
using Statistics
using Arrow 
using StructTypes
using JSON3
using DataFrames 
const PSID = PowerSimulationsDynamics
const PSY = PowerSystems

include("../src/train.jl")
include("../src/constants.jl")
include("../src/DynamicComponents.jl")
include("../src/init_functions.jl") #get rid of this? should need one function for all surrogates 
include("../src/instantiate.jl")
include("../src/SurrogateModels.jl")
include("../src/utils.jl")
include("../src/parameter_utils.jl")
include("../src/visualize.jl")
configure_logging(console_level = Logging.Error)


sys_train = System("systems/sys_train.json")
train_data = TrainData(sys_train, "data/train_input_data")
train_params = JSON3.read(read("data/default_NODE_params.json"), NODETrainParams)
train_params.export_mode = 1 

ll = train(train_params, train_data)
p = visualize_training(train_params, train_data)



## Example
# import("my_functions.jl")

# train_data = TrainData(sys_train, "test.arrow")
# train_params = UODETrainParams(Args[1])

# status = train(train_params, train_data)
################################################################ SLURM FILE 
# julia --project params_file.json