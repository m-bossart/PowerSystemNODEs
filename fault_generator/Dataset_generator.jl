using Pkg
Pkg.activate("Project.toml")
Pkg.instantiate
Pkg.status
using PowerSimulationsDynamics
PSID = PowerSimulationsDynamics
PSID.GeneratorTrip
using PowerSystems
using Sundials
using Plots
using YAML
using CSV
gr()

include("generator.jl")
include("fault.jl")

number = 1;
configuration = YAML.load_file("config.yml")
SimulationParameters = configuration["SimulationParameters"]
FaultParameters = configuration["FaultParameters"]
OutputParameters = configuration["OutputParameters"]
fault_types = configuration["FaultsType"]

system = System(configuration["CompleteSystemPath"])
time_span = (SimulationParameters["TspanStart"], SimulationParameters["TspanEnd"])

for fault in FaultParameters
    
    temp = get_fault(string(fault.first), FaultParameters, system)
    sim = PSID.Simulation(PSID.ResidualModel, system, pwd(), time_span, temp)
    PSID.execute!(
    sim, 
    IDA(),
    dtmax = SimulationParameters["DtMax"],
    );

    results = read_results(sim)
    angles = []
    volts = []
    for output in OutputParameters["OutputData"]
        append!(angles, get_state_series(results, (output, :δ));)
        append!(volts, get_voltage_magnitude_series(results, OutputParameters["BusNumber"]));
    end
    CSV.write(OutputParameters["OutputFile"] * string(fault.first) * ".csv", (x=volts, y=angles))

end