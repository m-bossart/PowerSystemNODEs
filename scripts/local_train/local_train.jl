using Revise
using PowerSimulationNODE
import PowerSimulationsDynamicsSurrogates
const PSIDS = PowerSimulationsDynamicsSurrogates
using Logging

#STUDY TIMING: 
######################################################################################
################################### SET PARAMETERS ###################################
######################################################################################
training_directory = "train_local"
p = TrainParams(
    base_path = training_directory,
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
            solver = "IDA",
            solver_tols = (1e-9, 1e-6),
            tspan = (0.0, 10.0),
            steps = 100,
            tsteps_spacing = "linear",
            formulation = "Residual",
            all_lines_dynamic = false,
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
            solver = "IDA",
            solver_tols = (1e-9, 1e-6),
            tspan = (0.0, 10.0),
            steps = 100,
            tsteps_spacing = "linear",
            formulation = "Residual",
            all_lines_dynamic = false,
            seed = 2,
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
            solver = "IDA",
            solver_tols = (1e-9, 1e-6),
            tspan = (0.0, 10.0),
            steps = 100,
            tsteps_spacing = "linear",
            formulation = "Residual",
            all_lines_dynamic = false,
            seed = 3,
        ),
    ),
    hidden_states = 10,
    model_initializer = (
        type = "dense",
        n_layer = 2,
        width_layers = 10,
        activation = "hardtanh",
    ),
    model_node = (
        type = "dense",
        n_layer = 1,
        width_layers = 10,
        activation = "hardtanh",
        σ2_initialization = 0.0,
    ),
    model_observation = (
        type = "dense",
        n_layer = 1,
        width_layers = 10,
        activation = "hardtanh",
    ),
    scaling_limits = (input_limits = (-1.0, 1.0), target_limits = (-1.0, 1.0)),
    steady_state_solver = (
        solver = "SSRootfind",  ##SSRootfind, Rodas4
        abstol = 1e-6, #1e-4        #xtol, ftol  #High tolerance -> standard NODE with initializer and observation 
        maxiters = 50,  #does not currently work for NLsolve (not set appropriately, always defaults to 1000 )
    ),
    optimizer = (
        sensealg = "Zygote",
        primary = "Adam",
        primary_η = 0.00000000000001,
        adjust = "nothing",
        adjust_η = 0.0,
    ),
    dynamic_solver = (solver = "Rodas5", tols = (1e-6, 1e-3), maxiters = 1e3), #1e-9, 1e-6 
    maxiters = 5, #TODO modify 
    lb_loss = 0.0,
    curriculum = "none",
    curriculum_timespans = [(tspan = (0.0, 10.0), batching_sample_factor = 1.0)],
    validation_loss_every_n = 100, #TODO modify 
    loss_function = (
        component_weights = (
            initialization_weight = 1.0,
            dynamic_weight = 1.0,
            residual_penalty = 1.0e9,
        ),
        type_weights = (rmse = 1.0, mae = 0.0),
    ),
    rng_seed = 1,
    output_mode_skip = 1,
    train_time_limit_seconds = 1e9,
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
##
animate_training(input_param_file, skip = 1)
a = generate_summary(joinpath(p.base_path, "output_data"))
pp = visualize_summary(a)
print_high_level_output_overview(a, p.base_path)

#= sys = System(joinpath(pwd(), "systems", "IEEE_14bus_modified.json"))
sys_train = System(joinpath(pwd(), "train_11", "system_data", "train_system.json"))
sys_validation = System(joinpath(pwd(), "train_11", "system_data", "validation_system.json")) =#
