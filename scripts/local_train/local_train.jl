using Revise
using PowerSimulationNODE
import PowerSimulationsDynamicsSurrogates
const PSIDS = PowerSimulationsDynamicsSurrogates
using Logging
using Serialization
using Plots
include("../build_datasets/utils.jl")
train_folder = "train_local3"
system_name = "CTESN_18bus_modified"
project_folder = "PowerSystemNODEs"
scratch_path = joinpath(pwd(), "..")
#Copy the full system over to the training directory.
mkpath(
    joinpath(
        scratch_path,
        project_folder,
        train_folder,
        PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
    ),
)
cp(
    joinpath(scratch_path, project_folder, "systems", string(system_name, ".json")),
    joinpath(
        scratch_path,
        project_folder,
        train_folder,
        PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
        string(system_name, ".json"),
    ),
    force = true,
)
cp(
    joinpath(
        scratch_path,
        project_folder,
        "systems",
        string(system_name, "_validation_descriptors.json"),
    ),
    joinpath(
        scratch_path,
        project_folder,
        train_folder,
        PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
        string(system_name, "_validation_descriptors.json"),
    ),
    force = true,
)

#p = params_data[2]
#PowerSimulationNODE._rebase_path!(p, joinpath(scratch_path, project_folder, train_folder))

######################################################################################
################################### SET PARAMETERS ###################################
######################################################################################

p = TrainParams(
    base_path = joinpath(scratch_path, project_folder, train_folder),
    system_path = joinpath("systems", "CTESN_18bus_modified.json"),
    surrogate_buses = [20],
    train_data = (
        id = "1",
        operating_points = PSIDS.SurrogateOperatingPoint[PSIDS.GenerationLoadScale(
            generation_scale = 1.0,
            load_scale = 1.0,
        ),],
        perturbations = repeat(
            [[PSIDS.RandomLoadChange(time = 1.0, load_multiplier_range = (0.0, 2.0))]],
            5,
        ),
        params = PSIDS.GenerateDataParams(
            solver = "Rodas5",
            solver_tols = (reltol = 1e-3, abstol = 1e-6),
            tspan = (0.0, 10.0),
            tstops = 0.0:0.1:10.0,
            tsave = 0.0:0.1:10.0,
            formulation = "MassMatrix",
            all_branches_dynamic = false,
            all_lines_dynamic = true,
            seed = 1,
        ),
        system = "full",
    ),
    validation_data = (
        id = "1",
        operating_points = PSIDS.SurrogateOperatingPoint[PSIDS.GenerationLoadScale(
            generation_scale = 1.0,
            load_scale = 1.0,
        ),],
        perturbations = repeat(
            [[PSIDS.RandomLoadChange(time = 1.0, load_multiplier_range = (0.0, 2.0))]],
            5,
        ),
        params = PSIDS.GenerateDataParams(
            solver = "Rodas5",
            solver_tols = (reltol = 1e-3, abstol = 1e-6),
            tspan = (0.0, 10.0),
            tstops = 0.0:0.1:10.0,
            tsave = 0.0:0.1:10.0,
            formulation = "MassMatrix",
            all_branches_dynamic = false,
            all_lines_dynamic = true,
            seed = 1,
        ),
    ),
    test_data = (
        id = "1",
        operating_points = PSIDS.SurrogateOperatingPoint[PSIDS.GenerationLoadScale(
            generation_scale = 1.0,
            load_scale = 1.0,
        ),],
        perturbations = repeat(
            [[PSIDS.RandomLoadChange(time = 1.0, load_multiplier_range = (0.0, 2.0))]],
            5,
        ),
        params = PSIDS.GenerateDataParams(
            solver = "Rodas5",
            solver_tols = (reltol = 1e-3, abstol = 1e-6),
            tspan = (0.0, 10.0),
            tstops = 0.0:0.1:10.0,
            tsave = 0.0:0.1:10.0,
            formulation = "MassMatrix",
            all_branches_dynamic = false,       #possible with current version of PSID? 
            all_lines_dynamic = true,
            seed = 1,
        ),
    ),
    model_params = SteadyStateNODEObsParams(
        name = "source_1",
        n_ports = 1,
        initializer_layer_type = "dense",
        initializer_n_layer = 2,
        initializer_width_layers = 10,
        initializer_activation = "hardtanh",
        dynamic_layer_type = "dense",
        dynamic_hidden_states = 10,
        dynamic_n_layer = 2,
        dynamic_width_layers = 10,
        dynamic_activation = "hardtanh",
        dynamic_σ2_initialization = 0.0,
        observation_layer_type = "dense",
        observation_n_layer = 1,
        observation_width_layers = 10,
        observation_activation = "hardtanh",
    ),
    steady_state_solver = (
        solver = "SSRootfind",
        abstol = 1e-4,
        maxiters = 1,  #TODO does not currently work for NLsolve (not set appropriately, always defaults to 1000 )
    ),
    dynamic_solver = (
        solver = "Rodas5",
        reltol = 1e-3,
        abstol = 1e-6,
        maxiters = 1e5,
        force_tstops = true,
    ),
    optimizer = [
        (  #PRIMARY! 
            sensealg = "Zygote",
            algorithm = "Adam",
            η = 0.000000000001,
            initial_stepnorm = 0.0,
            maxiters = 5,
            lb_loss = 0.0,
            curriculum = "individual faults",
            curriculum_timespans = [(tspan = (0.0, 10.0), batching_sample_factor = 1.0)],
            fix_params = [:initializer],
            loss_function = (
                component_weights = (
                    initialization_weight = 1.0,
                    dynamic_weight = 1.0,
                    residual_penalty = 1.0e9,
                ),
                type_weights = (rmse = 1.0, mae = 0.0),
            ),
        ),
        (  #Secondary 
            sensealg = "Zygote",
            algorithm = "Bfgs",
            η = 0.0,
            initial_stepnorm = 0.001,
            maxiters = 5,
            lb_loss = 0.0,
            curriculum = "individual faults",
            curriculum_timespans = [(tspan = (0.0, 10.0), batching_sample_factor = 1.0)],
            fix_params = [],
            loss_function = (
                component_weights = (
                    initialization_weight = 1.0,
                    dynamic_weight = 1.0,
                    residual_penalty = 1.0e9,
                ),
                type_weights = (rmse = 1.0, mae = 0.0),
            ),
        ),
    ],
    p_start = [],
    validation_loss_every_n = 100, #TODO modify 
    rng_seed = 1,
    output_mode_skip = 1,
    train_time_limit_seconds = 1e9,
)
##
######################################################################################
################################# BUILD AND GENERATE #################################
######################################################################################
build_subsystems(p)
mkpath(joinpath(p.base_path, PowerSimulationNODE.INPUT_FOLDER_NAME))
generate_train_data(p)
generate_validation_data(p)
generate_test_data(p)

##########################
# Visualize datasets (use before attempting to train)
train_dataset = Serialization.deserialize(p.train_data_path)
visualize_dataset(train_dataset)
validation_dataset = Serialization.deserialize(p.validation_data_path)
visualize_dataset(validation_dataset)
test_dataset = Serialization.deserialize(p.test_data_path)
visualize_dataset(test_dataset)

######################################################################################
####################################### TRAIN ########################################
######################################################################################
train(p)
######################################################################################
############################### ANALYZE AND VISUALIZE ################################
######################################################################################
##
input_param_file = joinpath(p.base_path, "input_data", "input_test1.json")
PowerSimulationNODE.serialize(p, input_param_file)
visualize_training(input_param_file, skip = 1)
##
animate_training(input_param_file, skip = 1)
a = generate_summary(joinpath(p.base_path, "output_data"))
pp = visualize_summary(a)
print_high_level_output_overview(a, p.base_path)

#= sys = System(joinpath(pwd(), "systems", "IEEE_14bus_modified.json"))
sys_train = System(joinpath(pwd(), "train_11", "system_data", "train_system.json"))
sys_validation = System(joinpath(pwd(), "train_11", "system_data", "validation_system.json")) =#
