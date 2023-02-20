using PowerSimulationNODE
using PowerSimulationsDynamicsSurrogates
const PSIDS = PowerSimulationsDynamicsSurrogates
include(joinpath(@__DIR__, "utils.jl"))
if Sys.iswindows() || Sys.isapple()
    const SCRATCH_PATH = joinpath(pwd(), "..")
else
    const SCRATCH_PATH = "/scratch/alpine/mabo4366"
end
train_folder = "exp_data_random"    #The name of the folder where everything related to the group of trainings will be stored (inputs, outputs, systems, logging, etc.)
system_name = "36Bus"           #The specific system from the "systems" folder to use. Will be copied over to the train_folder to make it self-contained.
project_folder = "PowerSystemNODEs"

_copy_full_system_to_train_directory(
    SCRATCH_PATH,
    project_folder,
    train_folder,
    system_name,
)

base_option = TrainParams(
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
            15,
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
            15,
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
            15,
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
    model_params = SteadyStateNODEParams(
        name = "source_1",
        n_ports = 1,
        initializer_layer_type = "dense",
        initializer_n_layer = 2,
        initializer_width_layers = 10,
        initializer_activation = "hardtanh",
        dynamic_layer_type = "dense",
        dynamic_hidden_states = 5,
        dynamic_n_layer = 1,
        dynamic_width_layers = 10,
        dynamic_activation = "hardtanh",
        dynamic_σ2_initialization = 0.0,
    ),
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
            sensealg = "Zygote",
            algorithm = "Adam",
            η = 0.01,
            initial_stepnorm = 0.0,
            maxiters = 6000,
            lb_loss = 0.0,
            curriculum = "individual faults",
            curriculum_timespans = [(tspan = (0.0, 10.0), batching_sample_factor = 1.0)],
            fix_params = [],
            loss_function = (α = 0.5, β = 0.5, residual_penalty = 1.0e9),
        ),
    ],
    validation_loss_every_n = 50,
    rng_seed = 1,
    output_mode_skip = 1,
    train_time_limit_seconds = 1e9,
    base_path = joinpath(SCRATCH_PATH, project_folder, train_folder),
    system_path = joinpath(
        SCRATCH_PATH,
        project_folder,
        train_folder,
        PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
        string(system_name, ".json"),
    ),
)

total_runs = 10
r1 = (:rng_seed, (min = 1, max = 1000))
#MODEL PARAMTERS
r2 = (:initializer_n_layer, (min = 1, max = 5))
r3 = (:initializer_width_layers, (min = 5, max = 20))
r4 = (:initializer_activation, (min = "na", max = "na", set = ["relu"]))
r5 = (:dynamic_hidden_states, (min = 5, max = 30))
r6 = (:dynamic_n_layer, (min = 1, max = 5))
r7 = (:dynamic_width_layers, (min = 5, max = 20))
r8 = (:dynamic_activation, (min = "na", max = "na", set = ["relu"]))
r9 = (:dynamic_σ2_initialization, (min = "na", max = "na", set = [0.0]))
#OPTIMIZER PARAMETERS
r10 = (:η, (min = 0.0, max = 0.1))
r11 = (:α, (min = 0.0, max = 1.0))
r12 = (:β, (min = 0.0, max = 1.0))

params_data = build_random_search!(
    base_option,
    total_runs,
    r1,
    r2,
    r3,
    r4,
    r5,
    r6,
    r7,
    r8,
    r9,
    r10,
    r11,
    r12,
)

#=
 hpc_params = SavioHPCTrain(;
    username = "jdlara",
    params_data = params_data,
    project_folder = "PowerSystemNODEs",
    scratch_path = "/global/home/users/jdlara",
)
  =#
hpc_params = AlpineHPCTrain(;
    username = "mabo4366",
    params_data = params_data,
    project_folder = project_folder,
    train_folder = train_folder,
    scratch_path = SCRATCH_PATH,
    time_limit_train = "23:59:59",             #Options: ["00:30:00", "23:59:59"]
    time_limit_generate_data = "02:00:00",
    QoS = "normal",
    partition = "amilan",
    force_generate_inputs = true,
    mb_per_cpu = 9600,  #Avoide OOM error on HPC 
)
generate_train_files(hpc_params)
##                                   
run_parallel_train(hpc_params)
