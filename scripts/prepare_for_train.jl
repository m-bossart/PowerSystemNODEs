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
using StructTypes
using JSON3
using DataFrames
using Random
const PSID = PowerSimulationsDynamics
const PSY = PowerSystems

include("../system_data/dynamic_components_data.jl")
include("../src/PowerSystemNODEs.jl")
configure_logging(console_level = Logging.Error)

base_system_path = joinpath(INPUT_SYSTEM_FOLDER_NAME, "base_system_3invs.json")
sys_full = System(base_system_path)            #Has all dynamic models, ready to simulate 

#pvs_data = fault_data_generator(sys_full)      #TODO - generate pvs data from simulation
#pvs_data = XXX()                               #TODO - generate pvs data from system id 
#sys_pvs = build_pvs(pvs_data, smooth=true)     #TODO - build pvs from data 
#include("../scripts/generate_fault_profiles.jl")    #Replace generate_fault_profiles.jl with three lines above 

sys_pvs =
    System(joinpath(INPUT_SYSTEM_FOLDER_NAME, "fault_library_3invs_vsms_20%lossP.json"))
#ground_truth_data = fault_data_generator(sys_full) #TODO 

label_area!(sys_full, [16], "surrogate")
@assert check_single_connecting_line_condition(sys_full)
sys_surr = remove_area(sys_full, "1")
sys_train = build_train_system(sys_surr, sys_pvs, "surrogate")
to_json(sys_train, joinpath(INPUT_FOLDER_NAME, "sys_train.json"), force = true)

#CHECK: Should now be able to delete generate_input_data.jl 
#create and serialize train data using default parameters 
d = generate_train_data(sys_train, NODETrainDataParams())
serialize(d, joinpath(INPUT_FOLDER_NAME, "data.json"))

#serialize default train params 
serialize(NODETrainParams(), joinpath(INPUT_FOLDER_NAME, "train_instance_1.json"))

######### POST TRAIN GENERATE PREDICTION DATA ########
#= sys_rest = remove_area(sys_full, "surrogate")
sys_reduced = build_reduced_system(sys_rest, NODE, "1")
prediction_data = fault_data_generator(sys_reduced) 
final_loss = loss(ground_truth_data, prediction_data)=#
