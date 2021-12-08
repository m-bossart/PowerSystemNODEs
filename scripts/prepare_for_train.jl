#This script includes and example of all the system/parameter/file generation needed before a training run.
#See figure in Readme for visualization of workflow. 
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
using YAML
using Sundials
using StructTypes
using JSON3
using DataFrames
using Random
using FFTW
const PSID = PowerSimulationsDynamics
const PSY = PowerSystems

include("../system_data/dynamic_components_data.jl")
include("../src/PowerSystemNODEs.jl")
configure_logging(console_level = Logging.Info)



base_system_path = joinpath(INPUT_SYSTEM_FOLDER_NAME, "14bus_3invs.json")
pvs_system_path = joinpath(INPUT_SYSTEM_FOLDER_NAME, "14bus_3invs_pvs.json")

include("build_full_system.jl")   #Build base system with all dynamic models 
sys_full = System(base_system_path)
pvs_data = fault_data_generator("scripts/config.yml") #pvs_data could be synthetic or real data 
sys_pvs = build_pvs(pvs_data)
to_json(sys_pvs, pvs_system_path, force = true)
sys_pvs = System(pvs_system_path)
#ground_truth_data = fault_data_generator(sys_full) 

label_area!(sys_full, [16], "surrogate")
@assert check_single_connecting_line_condition(sys_full)
sys_surr = remove_area(sys_full, "1")
sys_train = build_train_system(sys_surr, sys_pvs, "surrogate")
to_json(sys_train, joinpath(INPUT_FOLDER_NAME, "system.json"), force = true)

##
#create and serialize train data using default parameters 
d = generate_train_data(sys_train, NODETrainDataParams())   #BUG - only works for single fault

#plot(d.data[:tsteps], d.data[:ir_ground_truth])
serialize(d, joinpath(INPUT_FOLDER_NAME, "data.json"))

#serialize default train params 
serialize(NODETrainParams(), joinpath(INPUT_FOLDER_NAME, "train_1.json"))

######### POST TRAIN GENERATE PREDICTION DATA ########
#= sys_rest = remove_area(sys_full, "surrogate")
sys_reduced = build_reduced_system(sys_rest, NODE, "1")
prediction_data = fault_data_generator(sys_reduced) 
final_loss = loss(ground_truth_data, prediction_data)=#
