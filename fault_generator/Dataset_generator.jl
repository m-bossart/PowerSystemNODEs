using Pkg
Pkg.activate("SIIPExamples.jl/Project.toml")
Pkg.instantiate
using PowerSimulationsDynamics
PSID = PowerSimulationsDynamics
using PowerSystems
using Sundials
using Plots
using YAML
gr()

include("generator.jl")
include("fault.jl")
file_dir = joinpath(
    "SIIPExamples.jl",
    "script",
    "4_PowerSimulationsDynamics_examples",
    "Data",
)

configuration = YAML.load_file("config.yml")
SimulationParameters = configuration["SimulationParameters"]
FaultParameters = configuration["FaultParameters"]
OutputParameters = configuration["OutputParameters"]
fault_type = configuration["FaultType"]

system = System(joinpath(file_dir, "omib_sys.json"))
time_span = (SimulationParameters["TspanStart"], SimulationParameters["TspanEnd"])
fault = get_fault(fault_type, FaultParameters)

sim = PSID.Simulation(PSID.ResidualModel, system, pwd(), time_span, fault)

PSID.execute!(
    sim, #simulation structure
    IDA(), #Sundials DAE Solver
    dtmax = SimulationParameters["DtMax"],
);

results = read_results(sim)

angle = get_state_series(results, ("generator-102-1", :δ));
Plots.plot(angle, xlabel = "time", ylabel = "rotor angle [rad]", label = "rotor angle")
Plots.savefig("try.png")

volt = get_voltage_magnitude_series(results, 102);
Plots.plot(volt, xlabel = "time", ylabel = "Voltage [pu]", label = "V_2")
Plots.savefig("try_2.png")
