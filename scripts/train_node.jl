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

configure_logging(console_level = Logging.Info)
#configure_logging(;filename = "train_node.log")

include("../src/PowerSystemNODEs.jl")
include("../system_data/dynamic_components_data.jl")

script_train_params_file = "input_data/train_instance_1.json"

train_params_file = isempty(ARGS) ? script_train_params_file : ARGS[1]
train_params = NODETrainParams(train_params_file)
status = train(train_params)

if train_params.graphical_report
    #visualize all training runs to see what worked well
    p = visualize_summary(NODETrainParams().output_data_path)
    #plot(p)
    #visualize individual training runs of interest
    plots = visualize_training(train_params_1)
    plot(plots[1], plots[2], layout = (1, 2))
    plots = visualize_training(train_params_2)
    plot(plots[1], plots[2], layout = (1, 2))
end
