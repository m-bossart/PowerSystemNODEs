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
include("../src/init_functions.jl") #get rid of this? should need one function for all surrogates 
include("../src/instantiate.jl")
include("../src/SurrogateModels.jl")
include("../src/utils.jl")
include("../src/parameter_utils.jl")
include("../src/visualize.jl")
configure_logging(console_level = Logging.Error)

train_params = JSON3.read(read("train_parameters/train_instance_1.json"), NODETrainParams)

ll = @time train(train_params)

#visualize_summary(pwd())   #finds the outputs folder, picks up the (total time, total_iterations, final loss) for each instance. Plots. 

#p = visualize_training(train_params)    #only input is train_params, give an error if you can't find the folder where the results should be.    

## Example
# import("my_functions.jl")

# train_params = UODETrainParams(Args[1])

# status = train(train_params)
################################################################ SLURM FILE 
# julia --project params_file.json
