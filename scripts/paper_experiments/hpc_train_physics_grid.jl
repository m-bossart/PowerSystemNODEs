using PowerSystems
using PowerFlows
using PowerSimulationNODE
using PowerSimulationsDynamicsSurrogates
using Serialization
const PSIDS = PowerSimulationsDynamicsSurrogates
include(joinpath(@__DIR__, "..", "hpc_train", "utils.jl"))
if Sys.iswindows() || Sys.isapple()
    const SCRATCH_PATH = joinpath(pwd(), "..")
else
    const SCRATCH_PATH = "/scratch/alpine/mabo4366"
end
train_folder = "exp_physics_grid"    #The name of the folder where everything related to the group of trainings will be stored (inputs, outputs, systems, logging, etc.)
system_name = "36bus_fix"           #The specific system from the "systems" folder to use. Will be copied over to the train_folder to make it self-contained.
project_folder = "PowerSystemNODEs"

if Sys.iswindows() || Sys.isapple()
    include(joinpath(pwd(), "scripts", "hpc_train", "utils.jl"))
    sys = System(joinpath("systems", string(system_name, ".json")))
    surrogate_buses = vcat(21:29, 21:39)
    θ = determine_p_start(sys, surrogate_buses)
    Serialization.serialize(joinpath("starting_parameters", "p_start_physics"), θ)
end

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
    model_params = MultiDeviceParams(name = "source_1"),
    optimizer = [
        (
            auto_sensealg = "ForwardDiff",
            algorithm = "Adam",
            log_η = -7.0,
            initial_stepnorm = 0.0,
            maxiters = 2010,
            steadystate_solver = (
                solver = "NLSolveJL",
                reltol = 1e-4,
                abstol = 1e-4,
                termination = "RelSafeBest",
            ),
            dynamic_solver = (
                solver = "Rodas5",
                sensealg = "QuadratureAdjoint",
                reltol = 1e-3,
                abstol = 1e-6,
                maxiters = 1e5,
                force_tstops = true,
            ),
            lb_loss = 0.0,
            curriculum = "individual faults",
            curriculum_timespans = [(tspan = (0.0, 10.0), batching_sample_factor = 1.0)],
            fix_params = [   #fix zip to constant impedance 
                #:P_fraction_1,
                #:P_fraction_2,
                #:P_fraction_3,
                #:Q_fraction_1,
                :base_power_zip,
                :max_active_power_Z,
                :max_active_power_I,
                :max_active_power_P,
                :max_reactive_power_Z,
                :max_reactive_power_I,
                :max_reactive_power_P,
                :base_power_gfl,
                :rated_voltage_gfl,
                :rated_current_gfl,
                #:Kp_p,
                #:Ki_p,
                :ωz_gfl,
                #:Kp_q,
                #:Ki_q,
                :ωf_gfl,
                :kpc_gfl,
                :kic_gfl,
                :kffv_gfl,
                :voltage_gfl,
                :ω_lp,
                :kp_pll,
                :ki_pll,
                :lf_gfl,
                :rf_gfl,
                :cf_gfl,
                :lg_gfl,
                :rg_gfl,
                :base_power_gfm,
                :rated_voltage_gfm,
                :rated_current_gfm,
                #:Rp,    
                :ωz_gfm,
                #:kq,
                :ωf_gfm,
                :kpv,
                :kiv,
                :kffv_gfm,
                :rv,
                :lv,
                :kpc_gfm,
                :kic_gfm,
                :kffi,
                :ωad,
                :kad,
                :voltage_gfm,
                :lf_gfm,
                :rf_gfm,
                :cf_gfm,
                :lg_gfm,
                :rg_gfm,
            ],
            loss_function = (α = 0.5, β = 1.0, residual_penalty = 1.0e2),
        ),
    ],
    check_validation_loss_iterations = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000], #, 2000, 3000, 4000],
    p_start = Serialization.deserialize(joinpath("starting_parameters", "p_start_physics")),
    final_validation_loss = true,
    time_limit_buffer_seconds = 7200,
    rng_seed = 11,
    output_mode_skip = 100,
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


g1 = (:log_η, tuple(-7.0:0.2:-1.0 ...,))
params_data = build_grid_search!(base_option, g1);

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
    time_limit_train = "1-23:59:59",
    time_limit_generate_data = "0-02:00:00",
    QoS = "long",
    partition = "amilan",
    train_folder_for_data = nothing, #"data",
    mb_per_cpu = 9600,  #Avoide OOM error on HPC 
)
generate_train_files(hpc_params)
##                                   
run_parallel_train(hpc_params)
