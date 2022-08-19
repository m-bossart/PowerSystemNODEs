using Revise
using PowerSimulationNODE
import PowerSimulationsDynamicsSurrogates
const PSIDS = PowerSimulationsDynamicsSurrogates
using Logging

######################################################################################
################################### SET PARAMETERS ###################################
######################################################################################
training_directory = "train_100"
p = TrainParams(
    base_path = training_directory,
    system_path = joinpath("systems", "IEEE_14bus_modified.json"),
    surrogate_buses = [16],
    train_data = (
        operating_points = PSIDS.SurrogateOperatingPoint[PSIDS.GenerationLoadScale(
            generation_scale = 1.0,
            load_scale = 1.0,
        )],
        perturbations = [[
            PSIDS.RandomLoadChange(time = 0.5, load_multiplier_range = (0.9, 1.1)),
        ]],
        params = PSIDS.GenerateDataParams(
            solver = "Rodas4",
            formulation = "MassMatrix",
            solver_tols = (1e-4, 1e-4),
        ),
        system = "full",
    ),
    validation_data = (
        operating_points = PSIDS.SurrogateOperatingPoint[PSIDS.GenerationLoadScale(
            generation_scale = 1.0,
            load_scale = 1.0,
        )],
        perturbations = [[
            PSIDS.RandomLoadChange(time = 1.0, load_multiplier_range = (0.9, 1.1)),
        ]],
        params = PSIDS.GenerateDataParams(
            solver = "Rodas4",
            formulation = "MassMatrix",
            solver_tols = (1e-2, 1e-2),
        ),
    ),
    test_data = (
        operating_points = PSIDS.SurrogateOperatingPoint[PSIDS.GenerationLoadScale(
            generation_scale = 1.0,
            load_scale = 1.0,
        )],
        perturbations = [[
            PSIDS.RandomLoadChange(time = 0.5, load_multiplier_range = (0.9, 1.1)),
        ]],
        params = PSIDS.GenerateDataParams(
            solver = "Rodas4",
            formulation = "MassMatrix",
            solver_tols = (1e-4, 1e-4),
        ),
    ),
)

######################################################################################
################################# BUILD AND GENERATE #################################
######################################################################################
build_subsystems(p)
mkpath(joinpath(p.base_path, PowerSimulationNODE.INPUT_FOLDER_NAME))
generate_train_data(p)
generate_validation_data(p)
generate_test_data(p)

######################################################################################
####################################### TRAIN ########################################
######################################################################################
train(p)

######################################################################################
############################### ANALYZE AND VISUALIZE ################################
######################################################################################

input_param_file = joinpath(p.base_path, "input_data", "input_test1.json")
PowerSimulationNODE.serialize(p, input_param_file)
visualize_training(input_param_file, skip = 1)
animate_training(input_param_file, skip = 1)
a = generate_summary(joinpath(p.base_path, "output_data"))
pp = visualize_summary(a)
print_high_level_output_overview(a, p.base_path)

#= sys = System(joinpath(pwd(), "systems", "IEEE_14bus_modified.json"))
sys_train = System(joinpath(pwd(), "train_11", "system_data", "train_system.json"))
sys_validation = System(joinpath(pwd(), "train_11", "system_data", "validation_system.json")) =#
