using PowerSimulationNODE
using PowerSimulationsDynamicsSurrogates
const PSIDS = PowerSimulationsDynamicsSurrogates
include(joinpath(@__DIR__, "utils.jl"))
if Sys.iswindows() || Sys.isapple()
    const SCRATCH_PATH = joinpath(pwd(), "..")
else
    const SCRATCH_PATH = "/scratch/alpine/mabo4366"
end
train_folder = "exp_expository"    #The name of the folder where everything related to the group of trainings will be stored (inputs, outputs, systems, logging, etc.)
system_name = "CTESN_18bus_modified"           #The specific system from the "systems" folder to use. Will be copied over to the train_folder to make it self-contained.
project_folder = "PowerSystemNODEs"

_copy_full_system_to_train_directory(
    SCRATCH_PATH,
    project_folder,
    train_folder,
    system_name,
)

base_option = TrainParams(
    train_id = "BASE",
    surrogate_buses = [20],
    train_data = (
        id = "1",
        operating_points = PSIDS.SurrogateOperatingPoint[PSIDS.GenerationLoadScale(
            generation_scale = 1.0,
            load_scale = 1.0,
        ),],
        perturbations = repeat(
            [[PSIDS.RandomLoadChange(time = 1.0, load_multiplier_range = (0.0, 2.0))]],
            1,
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
            1,
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
            1,
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
        initializer_width_layers_relative_input = 10,
        initializer_activation = "tanh",
        dynamic_layer_type = "dense",
        dynamic_hidden_states = 5,
        dynamic_n_layer = 1,
        dynamic_width_layers_relative_input = 10,
        dynamic_activation = "tanh",
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
            log_η = -2.0,
            initial_stepnorm = 0.0,
            maxiters = 6000,
            lb_loss = 0.0,
            curriculum = "individual faults",
            curriculum_timespans = [(tspan = (0.0, 10.0), batching_sample_factor = 1.0)],
            fix_params = [],
            loss_function = (α = 0.5, β = 0.5, residual_penalty = 1.0e9),
        ),
    ],
    check_validation_loss_iterations = collect(1000:50:6000),
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

g1 = (:rng_seed, (min = 1, max = 1000))
#MODEL PARAMTERS
g2 = (:initializer_n_layer, (1, 5))
g3 = (:initializer_width_layers_relative_input, (5, 20))
g4 = (:dynamic_hidden_states, (5, 30))
g5 = (:dynamic_n_layer, (1, 5))
g6 = (:dynamic_width_layers_relative_input, (5, 20))
#OPTIMIZER PARAMETERS
g7 = (:log_η, (-5.0, -4.0, -3.0, -2.0))
g8 = (:α, (0.0, 0.33, 0.66, 1.0))
g9 = (:β, (0.0, 0.33, 0.66, 1.0))
params_data = build_grid_search!(base_option, g7, g8, g9)

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
    train_folder_for_data = nothing,
    mb_per_cpu = 9600,  #Avoide OOM error on HPC 
)
generate_train_files(hpc_params)
##                                   
run_parallel_train(hpc_params)
