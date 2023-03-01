#TODO - decide on base data-set.
using PowerSimulationNODE
using PowerSimulationsDynamicsSurrogates
const PSIDS = PowerSimulationsDynamicsSurrogates
include(joinpath(@__DIR__, "utils.jl"))
if Sys.iswindows() || Sys.isapple()
    const SCRATCH_PATH = joinpath(pwd(), "..")
else
    const SCRATCH_PATH = "/scratch/alpine/mabo4366"
end
train_folder = "exp_physics_random"    #The name of the folder where everything related to the group of trainings will be stored (inputs, outputs, systems, logging, etc.)
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
            log_η = -2.0,
            initial_stepnorm = 0.0,
            maxiters = 6000,
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
    p_start = Float64[
        0.2,
        0.2,
        0.4,
        0.4,
        0.4,
        0.4,
        1.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        12.700000000000001,
        615.9420289855071,
        1.9835456702225953,
        30.679163077003263,
        41.818017263110214,
        3.899221794528066,
        30.29304511089826,
        42.88608307322704,
        5.813096036884919,
        11.866710997978586,
        0.0,
        600.0,
        512.5954002113241,
        0.08625118210709415,
        4.482557140304551,
        0.07999999999999999,
        0.0029999999999999996,
        0.07399999999999998,
        0.19999999999999998,
        0.009999999999999998,
        12.700000000000001,
        205.3140096618357,
        0.04991162961714421,
        30.41129394073096,
        0.2018427749817273,
        1006.7417873074656,
        0.579584126157358,
        752.8680589602519,
        0.0,
        0.0,
        0.19950887723724778,
        1.2844250396816308,
        1.2844250396816308,
        0.0,
        47.58554332620469,
        0.1975149812011315,
        600.0,
        0.07999999999999999,
        0.0029999999999999996,
        0.07399999999999998,
        0.19999999999999998,
        0.009999999999999998,
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

total_runs = 10
r_1 = (:log_η, (min = -8.0, max = -2.0))
r_2 = (:β, (min = 0.0, max = 1.0))
params_data = build_random_search!(base_option, total_runs, r_1, r_2)

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
