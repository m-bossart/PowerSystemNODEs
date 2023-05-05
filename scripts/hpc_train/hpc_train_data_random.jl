using PowerSimulationNODE
using PowerSimulationsDynamicsSurrogates
using Random
const PSIDS = PowerSimulationsDynamicsSurrogates
include(joinpath(@__DIR__, "utils.jl"))
if Sys.iswindows() || Sys.isapple()
    const SCRATCH_PATH = joinpath(pwd(), "..")
else
    const SCRATCH_PATH = "/scratch/alpine/mabo4366"
end
train_folder = "exp_data_random"    #The name of the folder where everything related to the group of trainings will be stored (inputs, outputs, systems, logging, etc.)
system_name = "36bus_fix"               #The specific system from the "systems" folder to use. Will be copied over to the train_folder to make it self-contained.
project_folder = "PowerSystemNODEs"

_copy_full_system_to_train_directory(
    SCRATCH_PATH,
    project_folder,
    train_folder,
    system_name,
)

base_option = TrainParams(
    train_id = "BASE",
    surrogate_buses = vcat(21:29, 31:39),
    train_data = (
        id = "1",
        operating_points = repeat(
            [
                RandomOperatingPointXiao(
                    generator_voltage_range = (0.94, 1.06),
                    generator_power_range = (0.0, 1.0),
                    load_multiplier_range = (0.5, 1.5),
                ),
            ],
            23,
        ),
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
            seed = 11,
        ),
        system = "full",
    ),
    validation_data = (
        id = "1",
        operating_points = repeat(
            [
                RandomOperatingPointXiao(
                    generator_voltage_range = (0.94, 1.06),
                    generator_power_range = (0.0, 1.0),
                    load_multiplier_range = (0.5, 1.5),
                ),
            ],
            20,
        ),
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
            seed = 22,
        ),
    ),
    test_data = (
        id = "1",
        operating_points = repeat(
            [
                RandomOperatingPointXiao(
                    generator_voltage_range = (0.94, 1.06),
                    generator_power_range = (0.0, 1.0),
                    load_multiplier_range = (0.5, 1.5),
                ),
            ],
            100,
        ),
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
            seed = 33,
        ),
    ),
    model_params = SteadyStateNODEParams(
        name = "source_1",
        n_ports = 1,
        initializer_layer_type = "dense",
        initializer_n_layer = 3,
        initializer_width_layers_relative_input = 15,
        initializer_activation = "tanh",
        dynamic_layer_type = "dense",
        dynamic_hidden_states = 10,
        dynamic_n_layer = 3,
        dynamic_width_layers_relative_input = 15,
        dynamic_activation = "tanh",
        dynamic_σ2_initialization = 0.0,
    ),
    optimizer = [
        (
            sensealg = "Zygote",
            algorithm = "Adam",
            log_η = -7.0,
            initial_stepnorm = 0.0,
            maxiters = 1000,
            steadystate_solver = (solver = "SSRootfind", abstol = 1e-4),
            dynamic_solver = (
                solver = "Rodas5",
                reltol = 1e-3,
                abstol = 1e-6,
                maxiters = 1e5,
                force_tstops = true,
            ),
            lb_loss = 0.0,
            curriculum = "individual faults",
            curriculum_timespans = [(tspan = (0.0, 10.0), batching_sample_factor = 1.0)],
            fix_params = [],
            loss_function = (α = 0.5, β = 1.0, residual_penalty = 1.0e2),
        ),
    ],
    check_validation_loss_iterations = [],
    final_validation_loss = true, 
    time_limit_buffer_seconds = 7200,
    rng_seed = 11,
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

total_runs = 100
#r1 = (:rng_seed, (min = 1, max = 1000))
r1 = (:initializer_n_layer, (min = 1, max = 2))
r2 = (:initializer_width_layers_relative_input, (min = 0, max = 15))
r3 = (:dynamic_n_layer, (min = 1, max = 2))
r4 = (:dynamic_hidden_states, (min = 2, max = 10))
r5 = (:dynamic_width_layers_relative_input, (min = 0, max = 15))
r6 = (:log_η, (min = -3.0, max = -1.0))
r7 = (:α, (min = 0.1, max = 0.9))

Random.seed!(1234)
params_data = build_random_search!(base_option, total_runs, r1, r2, r3, r4, r5, r6, r7)
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
    time_limit_train = "0-23:59:59",
    time_limit_generate_data = "0-02:00:00",
    QoS = "normal",
    partition = "amilan",
    train_folder_for_data = nothing,
    mb_per_cpu = 11250,  #Avoide OOM error on HPC?
)
generate_train_files(hpc_params)
##                                   
run_parallel_train(hpc_params)
