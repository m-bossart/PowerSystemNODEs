using Pkg
Pkg.activate(".")
Pkg.instantiate

using PowerSimulationsDynamics
#PSID = PowerSimulationsDynamics
using DifferentialEquations
using PowerSystems
using Sundials
using Plots
using YAML
using CSV
using DataFrames
using Arrow
include("fault.jl") #get_fault function 



df = fault_data_generator("fault_generator/config.yml")




