using PowerSystems
using PowerSimulationsDynamics
using PowerSimulationNODE
using PowerSimulationsDynamicsSurrogates
const PSIDS = PowerSimulationsDynamicsSurrogates
const PSN = PowerSimulationNODE
using Serialization
using Plots

include("build_datasets/utils.jl")
include("hpc_train/utils.jl")
train_folder = "exp_starting_parameters"
system_name = "36bus_fix"
project_folder = "PowerSystemNODEs"
scratch_path = joinpath(pwd(), "..")

_copy_full_system_to_train_directory(
    scratch_path,
    project_folder,
    train_folder,
    system_name,
)

sys = System(joinpath("systems", string(system_name, ".json")))
surrogate_buses = [21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37, 38, 39]

p_start = determine_p_start(sys, surrogate_buses);
