using PowerSimulationNODE
using PowerSimulationsDynamicsSurrogates
const PSIDS = PowerSimulationsDynamicsSurrogates

train_folder = "train_1"    #The name of the folder where everything related to the group of trainings will be stored (inputs, outputs, systems, logging, etc.)
system_name = "full_system" #The specific system from the "systems" folder to use. Will be copied over to the train_folder to make it self-contained.

#Copy the full system over to the training directory.
mkpath(joinpath(pwd(), train_folder, PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME))
cp(
    joinpath(pwd(), "systems", string(system_name, ".json")),
    joinpath(
        pwd(),
        train_folder,
        PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
        "system.json",
    ),
)
cp(
    joinpath(pwd(), "systems", string(system_name, "_validation_descriptors.json")),
    joinpath(
        pwd(),
        train_folder,
        PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
        "system_validation_descriptors.json",
    ),
)

params_data = TrainParams[]
no_change_params = Dict{Symbol, Any}()
change_params = Dict{Symbol, Any}()

#INDICATE CONSTANT, NON-DEFAULT PARAMETERS (surrogate_buses and system_path CANNOT change)
no_change_params[:maxiters] = 100
no_change_params[:surrogate_buses] = [2]
no_change_params[:base_path] = joinpath(pwd(), train_folder)
no_change_params[:train_data] =
    (
        id = "1",
        operating_points = PSIDS.SurrogateOperatingPoint[PSIDS.GenerationLoadScale()],
        perturbations = [[PSIDS.PVS(source_name = "InfBus")]],
        params = PSIDS.GenerateDataParams(),
        system = "reduced",     #generate from the reduced system with sources to perturb or the full system
    ),

    #INDICATE PARAMETES TO ITERATE OVER COMBINATORIALLY 
    change_params[:optimizer] = [
        (
            sensealg = "Zygote",
            primary = "Adam",
            primary_η = 0.0001,
            adjust = "nothing",
            adjust_η = 0.0,
        ),
        (
            sensealg = "Zygote",
            primary = "Adam",
            primary_η = 0.0001,
            adjust = "nothing",
            adjust_η = 0.0,
        ),
    ]

change_params[:test_data] = [
    (
        id = "1",
        operating_points = PSIDS.SurrogateOperatingPoint[PSIDS.GenerationLoadScale()],
        perturbations = [[PSIDS.PVS(source_name = "InfBus")]],
        params = PSIDS.GenerateDataParams(),
    ),
    (
        id = "2",
        operating_points = PSIDS.SurrogateOperatingPoint[PSIDS.GenerationLoadScale()],
        perturbations = [[PSIDS.PVS(source_name = "InfBus")]],
        params = PSIDS.GenerateDataParams(),
    ),
]

build_params_list!(params_data, no_change_params, change_params)
@warn "Number of trainings:", length(params_data)

#=
 hpc_params = SavioHPCTrain(;
    username = "jdlara",
    params_data = params_data,
    project_folder = "PowerSystemNODEs",
    scratch_path = "/global/home/users/jdlara",
)
  =#

hpc_params = SummitHPCTrain(;
    username = "mabo4366",
    params_data = params_data,
    project_folder = "PowerSystemNODEs",
    train_folder = train_folder,
    scratch_path = "/scratch/summit/mabo4366", #Options: [pwd(), "/scratch/summit/mabo4366"]
    time_limit_train = "00:30:00",             #Options: ["00:30:00", "23:59:59"]
    time_limit_generate_data = "00:30:00",
    QoS = "normal",
    partition = "shas-testing",                #Options: ["shas-testing", "shas"]
    force_generate_inputs = true,
    mb_per_cpu = 4800,
)

generate_train_files(hpc_params)
##                                   
run_parallel_train(hpc_params)
