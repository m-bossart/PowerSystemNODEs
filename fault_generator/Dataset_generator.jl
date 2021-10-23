using Pkg
Pkg.activate("files/SIIPExamples.jl/Project.toml")
Pkg.instantiate
using PowerSimulationsDynamics
PSID = PowerSimulationsDynamics
using PowerSystems
using Sundials
using Plots
using YAML


include("fault.jl")

configuration = YAML.load_file("config.yml")
SimulationParameters = configuration["SimulationParameters"]
FaultParameters = configuration["FaultParameters"]
OutputParameters = configuration["OutputParameters"]
fault_type = configuration["FaultType"]

system = System(configuration["CompleteSystemPath"])

time_span = (SimulationParameters["TspanStart"], SimulationParameters["TspanEnd"])

fault = get_fault(fault_type, FaultParameters)

sim = PSID.simulation(PSID.ImplicitModel, system, pwd(), time_span, fault)


PSID.execute!(
    sim, #simulation structure
    IDA(), #Sundials DAE Solver
    dtmax = SimulationParameters["DtMax"],
);

angle = get_state_series(sim, ("generator-102-1", :δ));
Plots.plot(angle, xlabel = "time", ylabel = "rotor angle [rad]", label = "rotor angle")