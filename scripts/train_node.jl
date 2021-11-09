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
using Random
const PSID = PowerSimulationsDynamics
const PSY = PowerSystems

include("../src/train.jl")
include("../src/constants.jl")
include("../src/DynamicComponents.jl")
include("../src/instantiate.jl")
include("../src/SurrogateModels.jl")
include("../src/utils.jl")
include("../src/visualize.jl")
configure_logging(console_level = Logging.Info)
#configure_logging(;filename = "train_node.log")

train_params_1 = JSON3.read(read("train_parameters/train_instance_1.json"), NODETrainParams)
status = train(train_params_1)

train_params_2 = JSON3.read(read("train_parameters/train_instance_1.json"), NODETrainParams)
train_params_2.train_id = "train_instance_2"
train_params_2.loss_function_scale = "none"
status = train(train_params_2)

#visualize all training runs to see what worked well
p = visualize_summary(NODETrainParams().output_data_path)
plot(p)

#visualize individual training runs of interest 
plots = visualize_training(train_params_1)
plot(plots[1], plots[2], layout = (1, 2))
plots = visualize_training(train_params_2)
plot(plots[1], plots[2], layout = (1, 2))

## Example
# import("my_functions.jl")

# train_params = UODETrainParams(Args[1])

# status = train(train_params)
################################################################ SLURM FILE 
# julia --project params_file.json
