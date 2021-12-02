#using Pkg
#Pkg.activate("Project.toml")
#Pkg.instantiate
using PowerSimulationsDynamics
PSID = PowerSimulationsDynamics
using PowerSystems
using Sundials
using Plots
using YAML

function fault_data_generator(path_to_config)
    configuration = YAML.load_file(path_to_config)

    SimulationParameters = configuration["SimulationParameters"]
    FaultParameters = configuration["FaultParameters"]
    OutputParameters = configuration["OutputParameters"]
    
    solver = SimulationParameters["Solver"] # TODO : instantiate solver 
    abstol = SimulationParameters["AbsTol"]
    reltol = SimulationParameters["RelTol"]
    tspan = (SimulationParameters["TspanStart"], SimulationParameters["TspanEnd"])
    step_size =  SimulationParameters["StepSize"]
    t_fault = SimulationParameters["FaultTime"]
    tsteps = tspan[1]:step_size:tspan[2]
    system = System(configuration["CompleteSystemPath"])
    
    faults = [] 
    append_faults!(faults, FaultParameters, system, t_fault) #Build list of PSID faults based on FaultParameters
    df = build_fault_data_dataframe(faults, system, OutputParameters)  #Run sim for each fault and build dataframe based on OutputParameters
    
    if OutputParameters["WriteFile"] == true  
        #= open(OutputParameters["OutputFile"], "w") do io
            Arrow.write(io, df[1]["data"])
        end =#
    end 
    return df 
end 


function build_fault_data_dataframe(faults, system, OutputParameters)
    output = Dict{Int, Dict{String,Any}}()
    for (i,fault) in enumerate(faults)  
        sim = Simulation!(MassMatrixModel, system, pwd(), tspan, fault)
        @warn fault
        execute!(
            sim,
            Rodas4(),  #TODO use from yaml 
            abstol = abstol,
            reltol = reltol,
            reset_simulation = false,
            saveat = tsteps,
        )
        results = read_results(sim)
    
        data = collect(tsteps)
        column_names = ["t"] 

        if OutputParameters["OutputData"]["BusNumbers"] !== nothing 
            for (i,bus_number) in enumerate(OutputParameters["OutputData"]["BusNumbers"])
                if OutputParameters["OutputData"]["BusData"][i] == "Vm"
                    push!(column_names, string(  bus_number, "_Vm"))
                    data = hcat(data, get_voltage_magnitude_series(results, bus_number)[2])
                elseif OutputParameters["OutputData"]["BusData"][i] == "Vtheta"
                    push!(column_names, string(  bus_number, "_Vtheta"))
                    data = hcat(data, get_voltage_angle_series(results, bus_number)[2])
                else 
                    @error "Invalid bus data, must be Vm or Vtheta"
                end 
            end
        end 
        if OutputParameters["OutputData"]["DynamicDevices"] !== nothing 
            for (i,device_name) in enumerate(OutputParameters["OutputData"]["DynamicDevices"])
                state_symbol = Symbol(OutputParameters["OutputData"]["States"][i])
                push!(column_names, string(  device_name, "_", state_symbol))
                data = hcat(data, get_state_series(results, (device_name, state_symbol))[2])
            end
        end 
        d = Dict("fault" => fault, "data"=>DataFrame(data,Symbol.(column_names))) 
        output[i] = d  

    end  
    return output 
end 





function append_faults!(faults, FaultParameters, system, t_fault)
    for fault_type in FaultParameters
        if fault_type[2]["DeviceName"] !== nothing 
            if fault_type[2]["DeviceName"] == "all"
                all_names = get_all_names(fault_type, system)
                for fault_device_name in all_names
                    f = get_fault(fault_device_name, fault_type, system, t_fault)
                    push!(faults, f)
                end  
            else 
                for fault_device_name in fault_type[2]["DeviceName"]
                    f = get_fault(fault_device_name, fault_type, system, t_fault)
                    push!(faults, f)
                end 
            end
        end 
    end   
end 


function get_all_names(fault_type, system)
    if fault_type[1] == "BranchTrip"
        names = get_name.(collect(get_components(Line, system)))
    elseif fault_type[1] == "BranchImpedanceChange"
        names = get_name.(collect(get_components(Line, system)))
    elseif fault_type[1] == "ControlReferenceChange"
        names = get_name.(collect(get_components(DynamicInjection, system)))
    elseif fault_type[1] == "GeneratorTrip"
        names = get_name.(collect(get_components(DynamicInjection, system)))
    elseif fault_type[1] == "LoadChange"
        names = get_name.(collect(get_components(ElectricLoad, system)))
    elseif fault_type[1] == "LoadTrip"
        names = get_name.(collect(get_components(ElectricLoad, system)))
    else
        @error "This type of fault is not supported"
    end
    return names 
end 


function get_fault(fault_device_name, fault_type, system, t_fault)
    #fault_name = fault_parameters[fault_type]
    if fault_type[1] == "BranchTrip"
        fault = BranchTrip(t_fault, Line, fault_device_name)    #Needs to be line?

    elseif fault_type[1] == "BranchImpedanceChange"
        fault = BranchImpedanceChange(t_fault, Line, fault_device_name, fault_type[2]["Multiplier"]) 

    elseif fault_type[1] == "ControlReferenceChange"
        g = get_component(DynamicInjection, system, fault_device_name)
        #TODO Get current reference, multiple by RefValue
        fault = ControlReferenceChange(t_fault, g, :P_ref, fault_type[2]["RefValue"])

    elseif fault_type[1] == "GeneratorTrip"
        g = get_component(DynamicInjection, system, fault_device_name)
        fault = GeneratorTrip(t_fault, g)

    elseif fault_type[1] == "LoadChange"
        l = get_component(ElectricLoad, system, fault_device_name)
        #TODO Get current reference, multiple by RefValue
        fault = PowerSimulationsDynamics.LoadChange(t_fault, l, :P_ref,fault_type[2]["RefValue"])

    elseif fault_type[1] == "LoadTrip"
        l = get_component(ElectricLoad, system, fault_device_name)
        fault = LoadTrip(t_fault, l)

    else
        error("This type of fault is not supported")
    end
    return fault
end

