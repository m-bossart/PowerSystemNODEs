using Pkg
Pkg.activate("Project.toml")
Pkg.instantiate
using PowerSimulationsDynamics
PSID = PowerSimulationsDynamics
using PowerSystems
using Sundials
using Plots
using YAML

function get_fault(fault_type, fault_parameters)
    fault_name = fault_parameters[fault_type]
    if fault_type == "BranchTrip"
        fault = BranchTrip(fault_name["Time"], Line, fault_name["BranchName"],)
    elseif fault_type == "BranchImpedanceChange"
        fault = BranchImpedanceChange(fault_name["Time"], fault_name["BranchType"], fault_name["BranchName"], fault_name["Multiplier"])        
    elseif fault_type == "NetworkSwitch"
        fault = NetworkSwitch(fault_name["Time"], fault_name["YbusRectangular"])
    elseif fault_type == "ControlReferenceChange"
        fault = ControlReferenceChange(fault_name["Time"], fault_name["Device"], fault_name["Signal"], fault_name["RefValue"])
    elseif fault_type == "GeneratorTrip"
        fault = GeneratorTrip(fault_name["Time"], fault_name["Device"])
    elseif fault_type == "SourceBusVoltageChange"
        fault = SourceBusVoltageChange(fault_name["Time"], fault_name["Device"], fault_name["SignalIndex"], fault_name["RefValue"])
    elseif fault_type == "LoadChange"
        fault = LoadChange(fault_name["Time"], fault_name["Device"], fault_name["Signal"], fault_name["RefValue"])
    elseif fault_type == "LoadTrip"
        fault = LoadTrip(fault_name["Time"], fault_name["Device"])
    else
        error("This type of fault is not supported")
    end
    return fault
end