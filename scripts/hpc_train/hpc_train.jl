using PowerSimulationNODE
using PowerSimulationsDynamicsSurrogates
const PSIDS = PowerSimulationsDynamicsSurrogates
if Sys.iswindows() || Sys.isapple()
    const SCRATCH_PATH = joinpath(pwd(), "..")
else
    const SCRATCH_PATH = "/scratch/alpine/mabo4366"
end
train_folder = "exp_1"    #The name of the folder where everything related to the group of trainings will be stored (inputs, outputs, systems, logging, etc.)
system_name = "36Bus" #The specific system from the "systems" folder to use. Will be copied over to the train_folder to make it self-contained.
project_folder = "PowerSystemNODEs"

#Copy the full system over to the training directory.
mkpath(
    joinpath(
        SCRATCH_PATH,
        project_folder,
        train_folder,
        PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
    ),
)
cp(
    joinpath(SCRATCH_PATH, project_folder, "systems", string(system_name, ".json")),
    joinpath(
        SCRATCH_PATH,
        project_folder,
        train_folder,
        PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
        string(system_name, ".json"),
    ),
    force = true,
)
cp(
    joinpath(
        SCRATCH_PATH,
        project_folder,
        "systems",
        string(system_name, "_validation_descriptors.json"),
    ),
    joinpath(
        SCRATCH_PATH,
        project_folder,
        train_folder,
        PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
        string(system_name, "_validation_descriptors.json"),
    ),
    force = true,
)

params_data = TrainParams[]
no_change_params = Dict{Symbol, Any}()
change_params = Dict{Symbol, Any}()

#INDICATE CONSTANT, NON-DEFAULT PARAMETERS (surrogate_buses and system_path CANNOT change; train_id set automatically)
no_change_params[:surrogate_buses] =
    [21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37, 38, 39]
no_change_params[:train_data] = (
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
)
no_change_params[:validation_data] = (
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
)
no_change_params[:test_data] = (
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
)
#= no_change_params[:model_params] = SteadyStateNODEParams(
    name = "source_1",
    n_ports = 1,
    initializer_layer_type = "dense", #Correct 
    initializer_n_layer = 2, #Correct 
    initializer_width_layers = 10, #Correct 
    initializer_activation = "hardtanh", #Correct 
    dynamic_layer_type = "dense", #Correct 
    dynamic_hidden_states = 5,  #Correct 
    dynamic_n_layer = 1, #Correct 
    dynamic_width_layers = 10, #Correct 
    dynamic_activation = "hardtanh", #Correct 
    dynamic_σ2_initialization = 0.0, #Correct 
) =#

#no_change_params[:steady_state_solver] = (solver = "SSRootfind", abstol = 1e-4)
no_change_params[:dynamic_solver] =
    (solver = "Rodas5", reltol = 1e-3, abstol = 1e-6, maxiters = 1e5, force_tstops = true)

no_change_params[:validation_loss_every_n] = 100

no_change_params[:output_mode_skip] = 1
no_change_params[:train_time_limit_seconds] = 1e9
no_change_params[:base_path] = joinpath(SCRATCH_PATH, project_folder, train_folder)
no_change_params[:system_path] = joinpath(
    SCRATCH_PATH,
    project_folder,
    train_folder,
    PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
    string(system_name, ".json"),
)

change_params[:rng_seed] = [1, 2]
change_params[:model_params] = [
    SteadyStateNODEParams(
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
    MultiDeviceParams(name = "source_1"),
]
change_params[:steady_state_solver] =
    [(solver = "SSRootfind", abstol = 1e-4), (solver = "SSRootfind", abstol = 1e10)]
change_params[:optimizer] = [
    [
        (
            sensealg = "Zygote",
            algorithm = "Adam",
            η = 0.01,
            initial_stepnorm = 0.0,
            maxiters = 2000,
            lb_loss = 0.0,
            curriculum = "individual faults",
            curriculum_timespans = [(tspan = (0.0, 10.0), batching_sample_factor = 1.0)],
            fix_params = [],
            loss_function = (
                component_weights = (
                    initialization_weight = 1.0,
                    dynamic_weight = 1.0,
                    residual_penalty = 1.0,
                ),
                type_weights = (rmse = 1.0, mae = 0.0),
            ),
        ),
    ],
    [
        (
            sensealg = "Zygote",
            algorithm = "Adam",
            η = 0.01,
            initial_stepnorm = 0.0,
            maxiters = 2000,
            lb_loss = 0.0,
            curriculum = "individual faults",
            curriculum_timespans = [(tspan = (0.0, 10.0), batching_sample_factor = 1.0)],
            fix_params = [
                :P_fraction_1,  #how P/Q is distributed among devices in MultiDevice
                :Q_fraction_1,
                :P_fraction_2,
                :Q_fraction_2,
                :P_fraction_3,
                :Q_fraction_3,
                :kffv_gfl,      #binary parameters in the inverters 
                :kffv_gfm,
                :kffi,
            ],
            loss_function = (
                component_weights = (
                    initialization_weight = 1.0,
                    dynamic_weight = 1.0,
                    residual_penalty = 1.0,
                ),
                type_weights = (rmse = 1.0, mae = 0.0),
            ),
        ),
    ],
]

#INDICATE PARAMETES TO ITERATE OVER COMBINATORIALLY 
build_params_list!(params_data, no_change_params, change_params)
@warn "Number of trainings:", length(params_data)
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
    force_generate_inputs = true,
    mb_per_cpu = 9600,  #Avoide OOM error on HPC 
)

generate_train_files(hpc_params)
##                                   
run_parallel_train(hpc_params)
