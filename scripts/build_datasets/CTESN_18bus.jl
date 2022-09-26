using PowerSystems
using PowerSimulationsDynamics
using PowerSimulationNODE
using PowerSimulationsDynamicsSurrogates
using Serialization
using Plots
const PSY = PowerSystems
const PSIDS = PowerSimulationsDynamicsSurrogates

include(joinpath(pwd(), "scripts", "build_datasets", "utils.jl"))
include(joinpath(pwd(), "system_data", "dynamic_components_data.jl"))

sys_name = "CTESN_18bus" #"CTESN_18bus" "CTESN_default_14bus"
surrogate_type = "gfl" #  ["constant_impedance", "gfm", "gfl"]
connecting_bus_number = 6 #13
plot_title = "$sys_name, surrogate @ $connecting_bus_number, surrogate $surrogate_type"

original_sys_path = joinpath(pwd(), "systems", string(sys_name, ".json"))
new_sys_path = joinpath(pwd(), "systems", string(sys_name, "_modified.json"))
sys = System(original_sys_path)

set_penetration_level!(sys, 0.4, 0.2)
sync_cap = 0.0
gfl_cap = 0.0
gfm_cap = 0.0
for c in get_components(
    DynamicGenerator{AndersonFouadMachine, SingleMass, AVRTypeI, TGTypeI, PSSFixed},
    sys,
)
    sync_cap += get_base_power(c)
end
for c in get_components(
    DynamicInverter{
        AverageConverter,
        OuterControl{ActivePowerDroop, ReactivePowerDroop},
        VoltageModeControl,
        FixedDCSource,
        FixedFrequency,
        LCLFilter,
    },
    sys,
)
    gfm_cap += get_base_power(c)
end
for c in get_components(
    DynamicInverter{
        AverageConverter,
        OuterControl{ActivePowerPI, ReactivePowerPI},
        CurrentModeControl,
        FixedDCSource,
        KauraPLL,
        LCLFilter,
    },
    sys,
)
    gfl_cap += get_base_power(c)
end
sync_cap / (sync_cap + gfl_cap + gfm_cap)
gfl_cap / (sync_cap + gfl_cap + gfm_cap)
gfm_cap / (sync_cap + gfl_cap + gfm_cap)

#Modifications - constant impedance loads 

for l in get_components(PowerLoad, sys)
    PSY.set_model!(l, LoadModels.ConstantImpedance)
end

#Add a surrogate part to the system 
display(solve_powerflow(sys)["bus_results"])
if surrogate_type == "constant_impedance"
    new_bus_number = add_surrogate_A!(sys, connecting_bus_number)
elseif surrogate_type == "gfm"
    new_bus_number = add_surrogate_B!(sys, connecting_bus_number)
elseif surrogate_type == "gfl"
    new_bus_number = add_surrogate_C!(sys, connecting_bus_number)
else
    @error "invalid surrogate type"
end
display(solve_powerflow(sys)["bus_results"])

to_json(sys, new_sys_path, force = true)
p = TrainParams(
    base_path = "test_dir",
    system_path = new_sys_path,
    surrogate_buses = [new_bus_number],
    validation_data = (
        id = "1",
        operating_points = PSIDS.SurrogateOperatingPoint[PSIDS.GenerationLoadScale(
            generation_scale = 1.0,
            load_scale = 1.0,
        ),],
        perturbations = repeat(
            [[PSIDS.RandomLoadChange(time = 1.0, load_multiplier_range = (0.0, 2.0))]],
            100,
        ),
        params = PSIDS.GenerateDataParams(
            solver = "Rodas5",
            solver_tols = (1e-9, 1e-6),
            tspan = (0.0, 10.0),
            steps = 1000,
            tsteps_spacing = "linear",
            formulation = "MassMatrix",
            all_lines_dynamic = false,
            seed = 1,
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
            10,
        ),
        params = PSIDS.GenerateDataParams(
            solver = "Rodas5",  #IDA
            solver_tols = (1e-9, 1e-6),
            tspan = (0.0, 10.0),
            steps = 1000,
            tsteps_spacing = "linear",
            formulation = "MassMatrix",
            all_lines_dynamic = true,
            seed = 2,
        ),
    ),
)

build_subsystems(p)
mkpath(joinpath("test_dir", "input_data"))  #add this to some function??
generate_validation_data(p)
dataset_validation = Serialization.deserialize(p.validation_data_path)
display(plot(visualize_dataset(dataset_validation), plot_title = plot_title))

#= generate_test_data(p)
dataset_test = Serialization.deserialize(p.test_data_path)
display(plot(visualize_dataset(dataset_test), plot_title = plot_title)) =#
