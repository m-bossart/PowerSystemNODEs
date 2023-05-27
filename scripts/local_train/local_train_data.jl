using Revise
using PowerSimulationNODE
using PowerSystems
using PowerSimulationsDynamics
using PowerSimulationsDynamicsSurrogates
const PSIDS = PowerSimulationsDynamicsSurrogates
using Logging
using Serialization
using Plots
include("../build_datasets/utils.jl")
include("../hpc_train/utils.jl")
train_folder = "train_local_data"
system_name = "36bus_fix"
project_folder = "PowerSystemNODEs"
scratch_path = joinpath(pwd(), "..")

_copy_full_system_to_train_directory(
    scratch_path,
    project_folder,
    train_folder,
    system_name,
)

#= for b in get_components(Bus, sys)
    if get_bustype(b) == BusTypes.PV ||  get_bustype(b) == BusTypes.REF
    println("bus: ", get_name(b))
    gens_at_bus = get_components(x-> get_bus(x) == b, ThermalStandard, sys)
    for g in gens_at_bus
        println("connected generator: ", get_name(g))
    end 
    loads_at_bus =  get_components(x-> get_bus(x) == b, StandardLoad, sys)
    for l in loads_at_bus
        println("connected load: ", get_name(l))
    end 
    println("")
end 
end 
 =#
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
        operating_points = repeat(
            [
                RandomOperatingPointXiao(
                    generator_voltage_range = (0.94, 1.06),  #(0.94, 1.06)
                    generator_power_range = (0.0, 1.0),      #(0.0, 1.0)
                    load_multiplier_range = (0.5, 1.5),      #(0.5, 1.5)
                ),
            ],
            1,
        ),
        perturbations = repeat(
            [[PSIDS.RandomLoadChange(time = 1.0, load_multiplier_range = (0.0, 2.0))]],
            23,
        ),
        params = PSIDS.GenerateDataParams(
            solver = "Rodas5", #"Rodas5
            solver_tols = (reltol = 1e-3, abstol = 1e-6),
            tspan = (0.0, 10.0),
            tstops = 0.0:0.1:10.0,
            tsave = 0.0:0.1:10.0,
            formulation = "MassMatrix", #"MassMatrix"
            all_branches_dynamic = false,
            all_lines_dynamic = false,
            seed = 11, #11
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
            1,
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
        initializer_n_layer = 1,
        initializer_width_layers_relative_input = 0,
        initializer_activation = "tanh",
        dynamic_layer_type = "dense",
        dynamic_hidden_states = 5,
        dynamic_n_layer = 1,
        dynamic_width_layers_relative_input = 0,
        dynamic_activation = "tanh",
        dynamic_σ2_initialization = 0.0,
    ),
    optimizer = [
        (
            sensealg = "Zygote",
            algorithm = "Adam",
            log_η = -9.0,
            initial_stepnorm = 0.0,
            maxiters = 25,        #CHanged to proof of concept 
            steadystate_solver = (solver = "Tsit5", abstol = 1e-4),
            dynamic_solver = (
                solver = "Rodas5",
                reltol = 1e-3,
                abstol = 1e-6,
                maxiters = Int64(1e5),
                force_tstops = true,
            ),
            lb_loss = 0.0,
            curriculum = "individual faults",
            curriculum_timespans = [(tspan = (0.0, 10.0), batching_sample_factor = 1.0)],
            fix_params = Symbol[],
            loss_function = (α = 0.5, β = 1.0, residual_penalty = 1.0e2),
        ),
    ],
    check_validation_loss_iterations = collect(1000:50:6000),
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
##
#= function _modify_test_data(p_path, new_test_data_params)
    p_old = TrainParams(p_path)
    p_old.test_data = new_test_data_params
    #generate_test_data(p_old)
    PowerSimulationNODE.serialize(p_old, p_path)
end

_modify_test_data(
    "transfers/exp_04_16_23_physics_grid/input_data/train_007.json",
    p2.test_data,
)        #This will modify the test data for the entire experiment. 
_modify_test_data(
    "transfers/exp_04_16_23_data_grid_long/input_data/train_006.json",
    p2.test_data,
)
 =#
######################################################################################
################################# BUILD AND GENERATE #################################
######################################################################################
#p = TrainParams("transfers/exp_04_16_23_physics_grid/input_data/train_007.json")
build_subsystems(p)
mkpath(joinpath(p.base_path, PowerSimulationNODE.INPUT_FOLDER_NAME))
generate_train_data(p)
generate_validation_data(p)
generate_test_data(p)
data_collection_location_validation =
    Serialization.deserialize(p.data_collection_location_path)[2]
##########################
# Visualize datasets (use before attempting to train)
train_dataset = Serialization.deserialize(p.train_data_path)
validation_dataset = Serialization.deserialize(p.validation_data_path)
display(visualize_dataset(validation_dataset))
test_dataset = Serialization.deserialize(p.test_data_path)
display(visualize_dataset(test_dataset))
_ = PowerSimulationNODE.instantiate_surrogate_flux(p, p.model_params, train_dataset)    #Just to instantiate the surrogate to see how many parameters it has. 

##############################################################################
####################################### TRAIN ########################################
######################################################################################
train(p)
######################################################################################
############################### ANALYZE AND VISUALIZE ################################
######################################################################################
#= ##
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
 =#
