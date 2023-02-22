using Revise
using PowerSimulationNODE
using PowerSimulationsDynamicsSurrogates
const PSIDS = PowerSimulationsDynamicsSurrogates
using Logging
using Serialization
using Plots
include("../build_datasets/utils.jl")
include("../hpc_train/utils.jl")
train_folder = "train_local_physics"
system_name = "36Bus"
project_folder = "PowerSystemNODEs"
scratch_path = joinpath(pwd(), "..")

_copy_full_system_to_train_directory(
    scratch_path,
    project_folder,
    train_folder,
    system_name,
)

######################################################################################
################################### SET PARAMETERS ###################################
######################################################################################
p = TrainParams(
    train_id = "BASE",
    surrogate_buses = [
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        31,
        32,
        33,
        34,
        35,
        36,
        37,
        38,
        39,
    ],
    train_data = (
        id = "1",
        operating_points = PSIDS.SurrogateOperatingPoint[
            PSIDS.GenerationLoadScale(generation_scale = 1.0, load_scale = 1.0),
            PSIDS.GenerationLoadScale(generation_scale = 0.9, load_scale = 0.9),
            PSIDS.GenerationLoadScale(generation_scale = 1.1, load_scale = 1.1),
        ],
        perturbations = repeat(
            [[PSIDS.RandomLoadChange(time = 1.0, load_multiplier_range = (0.0, 2.0))]],
            3,
        ),
        params = PSIDS.GenerateDataParams(
            solver = "Rodas5",
            solver_tols = (reltol = 1e-3, abstol = 1e-6),
            tspan = (0.0, 10.0),
            tstops = 0.0:0.1:10.0,
            tsave = 0.0:0.1:10.0,
            formulation = "MassMatrix",
            all_branches_dynamic = false,
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
            3,
        ),
        params = PSIDS.GenerateDataParams(
            solver = "Rodas5",
            solver_tols = (reltol = 1e-3, abstol = 1e-6),
            tspan = (0.0, 10.0),
            tstops = 0.0:0.1:10.0,
            tsave = 0.0:0.1:10.0,
            formulation = "MassMatrix",
            all_branches_dynamic = false,
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
            3,
        ),
        params = PSIDS.GenerateDataParams(
            solver = "Rodas5",
            solver_tols = (reltol = 1e-3, abstol = 1e-6),
            tspan = (0.0, 10.0),
            tstops = 0.0:0.1:10.0,
            tsave = 0.0:0.1:10.0,
            formulation = "MassMatrix",
            all_branches_dynamic = false,       #possible with current version of PSID? 
            all_lines_dynamic = false,
            seed = 3,
        ),
    ),
    model_params = MultiDeviceParams(name = "source_1"),
    steady_state_solver = (solver = "SSRootfind", abstol = 1e-4),
    dynamic_solver = (
        solver = "Rodas5",
        reltol = 1e-3,
        abstol = 1e-6,
        maxiters = 1e5,
        force_tstops = true,
    ),
    optimizer = [
        (
            sensealg = "ForwardDiff",
            algorithm = "Adam",
            log_η = -9.0,
            initial_stepnorm = 0.0,
            maxiters = 10,
            lb_loss = 0.0,
            curriculum = "individual faults",
            curriculum_timespans = [(tspan = (0.0, 10.0), batching_sample_factor = 1.0)],
            fix_params = [
                :P_fraction_1,
                :Q_fraction_1,
                :P_fraction_2,
                :Q_fraction_2,
                :P_fraction_3,
                :Q_fraction_3,
                :kffv_gfl,
                :kffv_gfm,
                :kffi,
            ],
            loss_function = (α = 0.5, β = 0.5, residual_penalty = 1.0e9),
        ),
    ],
    validation_loss_every_n = 50,
    rng_seed = 1,
    output_mode_skip = 1,
    train_time_limit_seconds = 1e9,
    base_path = joinpath(scratch_path, project_folder, train_folder),
    system_path = joinpath(
        scratch_path,
        project_folder,
        train_folder,
        PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
        string(system_name, ".json"),
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

##########################
# Visualize datasets (use before attempting to train)
train_dataset = Serialization.deserialize(p.train_data_path)
display(visualize_dataset(train_dataset))
validation_dataset = Serialization.deserialize(p.validation_data_path)
display(visualize_dataset(validation_dataset))
test_dataset = Serialization.deserialize(p.test_data_path)
display(visualize_dataset(test_dataset))

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
