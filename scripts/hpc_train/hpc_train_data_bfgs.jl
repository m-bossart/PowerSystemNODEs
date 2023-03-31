using PowerSimulationNODE
using PowerSimulationsDynamicsSurrogates
using Serialization
using JSON3
const PSIDS = PowerSimulationsDynamicsSurrogates
include(joinpath(@__DIR__, "utils.jl"))
if Sys.iswindows() || Sys.isapple()
    const SCRATCH_PATH = joinpath(pwd(), "..")
else
    const SCRATCH_PATH = "/scratch/alpine/mabo4366"
end
train_folder = "exp_data_bfgs"    #The name of the folder where everything related to the group of trainings will be stored (inputs, outputs, systems, logging, etc.)
system_name = "36Bus_CR"               #The specific system from the "systems" folder to use. Will be copied over to the train_folder to make it self-contained.
project_folder = "PowerSystemNODEs"

if Sys.iswindows() || Sys.isapple()
    starting_file =
        joinpath("transfers", "exp_03_30_23_data_grid", "input_data", "train_028.json")
    θ = _get_parameters_from_prior_training(starting_file)
    Serialization.serialize(joinpath("starting_parameters", "starting_parameters"), θ)
    p = TrainParams(starting_file)
    PowerSimulationNODE.serialize(
        p,
        joinpath("starting_parameters", "starting_parameters_trainparams"),
    )
end


p = TrainParams(joinpath("starting_parameters", "starting_parameters_trainparams"))
_copy_full_system_to_train_directory(
    SCRATCH_PATH,
    project_folder,
    train_folder,
    system_name,
)

base_option = TrainParams(
    train_id = "BASE",
    surrogate_buses = p.surrogate_buses, #vcat(21:29, 31:39),
    train_data = p.train_data,
    validation_data = p.validation_data,
    test_data = p.test_data,
    model_params = p.model_params,
    steady_state_solver = p.steady_state_solver,
    dynamic_solver = p.dynamic_solver,
    optimizer = [
        (
            sensealg = "Zygote",
            algorithm = "Bfgs",
            log_η = -2.0,
            initial_stepnorm = 0.01,
            maxiters = 100,
            lb_loss = 0.0,
            curriculum = "simultaneous",    #must be simultaneous for BFGS? 
            curriculum_timespans = [(tspan = (0.0, 10.0), batching_sample_factor = 1.0)],
            fix_params = [],
            loss_function = (α = 0.5, β = 1.0, residual_penalty = 1.0e2),
        ),
    ],
    p_start = Serialization.deserialize(
        joinpath("starting_parameters", "starting_parameters"),
    ),
    check_validation_loss_iterations = p.check_validation_loss_iterations,
    validation_loss_termination = p.validation_loss_termination,
    rng_seed = p.rng_seed,
    output_mode_skip = p.output_mode_skip,
    train_time_limit_seconds = p.train_time_limit_seconds,
    base_path = joinpath(SCRATCH_PATH, project_folder, train_folder),
    system_path = joinpath(
        SCRATCH_PATH,
        project_folder,
        train_folder,
        PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
        string(system_name, ".json"),
    ),
)

g = (:initial_stepnorm, (0.0, 0.001, 0.01, 0.1))

params_data = build_grid_search!(base_option, g)
##
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
    train_folder_for_data = "xiao_loadchange_30_10_10",
    mb_per_cpu = 9600,  #Avoide OOM error on HPC 
)
generate_train_files(hpc_params)
##                                   
run_parallel_train(hpc_params)
