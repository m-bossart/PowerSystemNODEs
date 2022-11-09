using PowerSimulationNODE
using PowerSimulationsDynamicsSurrogates
const PSIDS = PowerSimulationsDynamicsSurrogates

train_folder = "exp_1"    #The name of the folder where everything related to the group of trainings will be stored (inputs, outputs, systems, logging, etc.)
system_name = "CTESN_18bus_modified" #The specific system from the "systems" folder to use. Will be copied over to the train_folder to make it self-contained.
project_folder = "PowerSystemNODEs"
scratch_path = "/scratch/alpine/mabo4366"  #Options: [ joinpath(pwd(), ".."), "/scratch/alpine/mabo4366"]

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

params_data = TrainParams[]
no_change_params = Dict{Symbol, Any}()
change_params = Dict{Symbol, Any}()

#INDICATE CONSTANT, NON-DEFAULT PARAMETERS (surrogate_buses and system_path CANNOT change; train_id set automatically)
no_change_params[:surrogate_buses] = [20]
no_change_params[:train_data] =
    (
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
    )
    no_change_params[:validation_data] =
        (
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
                seed = 2,
            ),
        )
        no_change_params[:test_data] =
            (
                id = "1",
                operating_points = PSIDS.SurrogateOperatingPoint[PSIDS.GenerationLoadScale(
                    generation_scale = 1.0,
                    load_scale = 1.0,
                ),],
                perturbations = repeat(
                    [[
                        PSIDS.RandomLoadChange(
                            time = 1.0,
                            load_multiplier_range = (0.0, 2.0),
                        ),
                    ]],
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
                    seed = 3,
                ),
            )
            #no_change_params[:hidden_states] = 10
            no_change_params[:model_initializer] =
                (type = "dense", n_layer = 2, width_layers = 10, activation = "hardtanh")
no_change_params[:model_node] = (
    type = "dense",
    n_layer = 1,
    width_layers = 10,
    activation = "hardtanh",
    σ2_initialization = 0.0,
)
no_change_params[:model_observation] =
    (type = "dense", n_layer = 1, width_layers = 10, activation = "hardtanh")
no_change_params[:scaling_limits] =
    (input_limits = (-1.0, 1.0), target_limits = (-1.0, 1.0))
no_change_params[:steady_state_solver] =
    (solver = "SSRootfind", abstol = 1e-4, maxiters = 5)
no_change_params[:dynamic_solver] =
    (solver = "Rodas5", reltol = 1e-3, abstol = 1e-6, maxiters = 1e5)
#= no_change_params[:optimizer] = (
    sensealg = "Zygote",
    primary = "Adam",
    primary_η = 0.0001,
    adjust = "nothing",
    adjust_η = 0.0,
) =#
no_change_params[:lb_loss] = 0.0
no_change_params[:curriculum] = "none"
no_change_params[:curriculum_timespans] =
    [(tspan = (0.0, 10.0), batching_sample_factor = 1.0)]
no_change_params[:validation_loss_every_n] = 100
#= no_change_params[:loss_function] = (
    component_weights = (
        initialization_weight = 1.0,
        dynamic_weight = 1.0,
        residual_penalty = 1.0,
    ),
    type_weights = (rmse = 1.0, mae = 0.0),
) =#
no_change_params[:rng_seed] = 123
no_change_params[:output_mode_skip] = 10
no_change_params[:train_time_limit_seconds] = 1e9
no_change_params[:base_path] = joinpath(scratch_path, project_folder, train_folder)
no_change_params[:system_path] = joinpath(
    scratch_path,
    project_folder,
    train_folder,
    PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
    string(system_name, ".json"),
)

#INDICATE PARAMETES TO ITERATE OVER COMBINATORIALLY 
change_params[:hidden_states] = [4, 8, 12]
change_params[:loss_function] = [
    (
        component_weights = (
            initialization_weight = 1.0,
            dynamic_weight = 1.0,
            residual_penalty = 1.0e9,
        ),
        type_weights = (rmse = 1.0, mae = 0.0),
    ),
    (
        component_weights = (
            initialization_weight = 1.0,
            dynamic_weight = 10.0,
            residual_penalty = 1.0e9,
        ),
        type_weights = (rmse = 1.0, mae = 0.0),
    ),
    (
        component_weights = (
            initialization_weight = 1.0,
            dynamic_weight = 100.0,
            residual_penalty = 1.0e9,
        ),
        type_weights = (rmse = 1.0, mae = 0.0),
    ),
]

change_params[:optimizer] = [
    (
        sensealg = "Zygote",
        primary = "Adam",
        primary_η = 0.001,
        primary_maxiters = 2000,
        adjust = "Bfgs",
        adjust_η = 0.0,
        adjust_maxiters = 50, 
    ),
    (
        sensealg = "Zygote",
        primary = "Adam",
        primary_η = 0.0001,
        primary_maxiters = 2000,
        adjust = "Bfgs",
        adjust_η = 0.0,
        adjust_maxiters = 50, 
    ),
    (
        sensealg = "Zygote",
        primary = "Adam",
        primary_η = 0.00001,
        primary_maxiters = 2000,
        adjust = "Bfgs",
        adjust_η = 0.0,
        adjust_maxiters = 50, 
    ),
]

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
    scratch_path = scratch_path,
    time_limit_train = "23:59:59",             #Options: ["00:30:00", "23:59:59"]
    time_limit_generate_data = "01:00:00",
    QoS = "normal",
    partition = "amilan",
    force_generate_inputs = true,
    mb_per_cpu = 9600,  #Avoide OOM error on HPC 
)

generate_train_files(hpc_params)
##                                   
run_parallel_train(hpc_params)
