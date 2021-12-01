#using Pkg
#Pkg.activate("Project.toml")
#Pkg.instantiate
using PowerSimulationsDynamics
PSID = PowerSimulationsDynamics
using PowerSystems
using Sundials
using Plots
using YAML

function get_fault(fault_type, fault_parameters, system)
    fault_name = fault_parameters[fault_type]
    if fault_type == "BranchTrip"
        fault = BranchTrip(fault_name["Time"], Line, fault_name["BranchName"],)

    elseif fault_type == "BranchImpedanceChange"
        fault = BranchImpedanceChange(fault_name["Time"], Line, fault_name["BranchName"], fault_name["Multiplier"])

    elseif fault_type == "ControlReferenceChange"
        g = get_component(DynamicGenerator, system, fault_name["Device"])
        fault = ControlReferenceChange(fault_name["Time"], g, :P_ref, fault_name["RefValue"])

    elseif fault_type == "GeneratorTrip"
        g = get_component(DynamicGenerator, system, fault_name["Device"])
        fault = GeneratorTrip(fault_name["Time"], g)

    elseif fault_type == "LoadChange"
        l_device = get_component(ElectricLoad, system, fault_name["Device"])
        fault = PowerSimulationsDynamics.LoadChange(fault_name["Time"], l_device, :P_ref, fault_name["RefValue"])

    elseif fault_type == "LoadTrip"
        l_device = get_component(ElectricLoad, system, fault_name["Device"])
        fault = LoadTrip(fault_name["Time"], l_device)

    else
        error("This type of fault is not supported")
    end
    return fault
end