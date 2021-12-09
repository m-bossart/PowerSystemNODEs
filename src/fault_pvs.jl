#using Pkg
#Pkg.activate("Project.toml")
#Pkg.instantiate
using PowerSimulationsDynamics
using PowerSystems
using Sundials
using Plots
using YAML

function fault_data_generator(path_to_config)
    configuration = YAML.load_file(path_to_config)

    SimulationParameters = configuration["SimulationParameters"]
    FaultParameters = configuration["FaultParameters"]
    OutputParameters = configuration["OutputParameters"]

    t_fault = SimulationParameters["FaultTime"]
    system = node_load_system(configuration["CompleteSystemPath"])

    faults = []
    append_faults!(faults, FaultParameters, system, t_fault) #Build list of PSID faults based on FaultParameters
    df = build_fault_data_dataframe(faults, system, OutputParameters, SimulationParameters)  #Run sim for each fault and build dataframe based on OutputParameters

    if OutputParameters["WriteFile"] == true
        #= open(OutputParameters["OutputFile"], "w") do io
            Arrow.write(io, df[1]["data"])
        end =#
    end
    return df
end

function build_pvs(dict_fault_data; pad_signal = true)
    pad_signal

    sys_pvs = node_load_system(100.0)
    add_component!(
        sys_pvs,
        Bus(
            number = 1,
            name = "1",
            bustype = BusTypes.REF,
            angle = 0.0,
            magnitude = 1.0,
            voltage_limits = (min = 0.9, max = 1.1),
            base_voltage = 345.0,
        ),
    )
    pvs_bus = collect(get_components(Bus, sys_pvs, x -> x.bustype == BusTypes.REF))[1]

    for (key, value) in dict_fault_data
        t = value["data"][:, 1]
        Vm = value["data"][:, value["pvs"]["Vm_index"]]
        Vθ = value["data"][:, value["pvs"]["Vθ_index"]]
        if pad_signal
            t_new, Vm_new = pad_tanh(t, Vm)
            t_new, Vθ_new = pad_tanh(t, Vθ)
            t = t_new
            Vm = Vm_new
            Vθ = Vθ_new
        end

        #Dft 
        N = length(t)
        fs = (N - 1) / (t[end] - t[1])
        freqs = fftfreq(N, fs)
        freqs_pos = freqs[freqs .>= 0] * (2 * pi)

        F_V = fft(Vm)
        F_V = F_V[freqs .>= 0]
        F_V = F_V / N
        F_V[2:end] = F_V[2:end] * 2
        internal_voltage_coefficients = [(-imag(f), real(f)) for f in F_V[2:end]]

        F_θ = fft(Vθ)
        F_θ = F_θ[freqs .>= 0]
        F_θ = F_θ / N
        F_θ[2:end] = F_θ[2:end] * 2
        internal_angle_coefficients = [(-imag(f), real(f)) for f in F_θ[2:end]]

        inf_source = Source(
            name = string("source", string(key)),
            active_power = value["pvs"]["P"],
            available = false,
            reactive_power = value["pvs"]["Q"],
            bus = pvs_bus,
            R_th = 0.0,
            X_th = 5e-6,
            internal_voltage = Vm[1],
            internal_angle = Vθ[1],
        )

        fault_source = PeriodicVariableSource(
            name = get_name(inf_source),
            R_th = get_R_th(inf_source),
            X_th = get_X_th(inf_source),
            internal_voltage_bias = real(F_V[1]),
            internal_voltage_frequencies = freqs_pos[2:end],
            internal_voltage_coefficients = internal_voltage_coefficients,
            internal_angle_bias = real(F_θ[1]),
            internal_angle_frequencies = freqs_pos[2:end],
            internal_angle_coefficients = internal_angle_coefficients,
        )
        add_component!(sys_pvs, inf_source)
        add_component!(sys_pvs, fault_source, inf_source)
    end
    return sys_pvs
end

function build_fault_data_dataframe(faults, system, OutputParameters, SimulationParameters)
    solver = SimulationParameters["Solver"] # TODO : instantiate solver 
    abstol = SimulationParameters["AbsTol"]
    reltol = SimulationParameters["RelTol"]
    tspan = (SimulationParameters["TspanStart"], SimulationParameters["TspanEnd"])
    step_size = SimulationParameters["StepSize"]
    t_fault = SimulationParameters["FaultTime"]
    tsteps = tspan[1]:step_size:tspan[2]

    output = Dict{Int, Dict{String, Any}}()
    for (i, fault) in enumerate(faults)
        sim = Simulation!(
            MassMatrixModel,
            system,
            pwd(),
            tspan,
            fault,
            console_level = PSID_CONSOLE_LEVEL,
            file_level = PSID_FILE_LEVEL,
        )
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
        timeseries_data = collect(tsteps)
        column_names = ["t"]
        pvs_data = Dict()
        if OutputParameters["PVS"]["SavePVSData"]
            df = solve_powerflow(system)["flow_results"]
            pvs_branch_name = OutputParameters["PVS"]["PVSBranchName"]
            pvs_bus_number = OutputParameters["PVS"]["PVSBus"]
            #Grab the df row with the branch specified 
            pvs_connecting_branch_row = df[df.line_name .== pvs_branch_name, :]
            bus_from = pvs_connecting_branch_row[!, :bus_from][1]
            bus_to = pvs_connecting_branch_row[!, :bus_to][1]
            if bus_from == pvs_bus_number
                pvs_data["P"] = pvs_connecting_branch_row[!, :P_from_to][1]
                pvs_data["Q"] = pvs_connecting_branch_row[!, :Q_from_to][1]
            elseif bus_to == pvs_bus_number
                pvs_data["P"] = pvs_connecting_branch_row[!, :P_to_from][1]
                pvs_data["Q"] = pvs_connecting_branch_row[!, :Q_to_from][1]
            else
                @warn "PVS branch and bus aren't connected"
            end
        end
        if OutputParameters["OutputData"]["BusNumbers"] !== nothing
            for (i, bus_number) in enumerate(OutputParameters["OutputData"]["BusNumbers"])
                if OutputParameters["OutputData"]["BusData"][i] == "Vm"
                    push!(column_names, string(bus_number, "_Vm"))
                    if bus_number == pvs_bus_number
                        pvs_data["Vm_index"] = i + 1
                    end
                    timeseries_data = hcat(
                        timeseries_data,
                        get_voltage_magnitude_series(results, bus_number)[2],
                    )
                elseif OutputParameters["OutputData"]["BusData"][i] == "Vtheta"
                    push!(column_names, string(bus_number, "_Vtheta"))
                    if bus_number == pvs_bus_number
                        pvs_data["Vθ_index"] = i + 1
                    end
                    timeseries_data = hcat(
                        timeseries_data,
                        get_voltage_angle_series(results, bus_number)[2],
                    )
                else
                    @error "Invalid bus data, must be Vm or Vtheta"
                end
            end
        end
        if OutputParameters["OutputData"]["DynamicDevices"] !== nothing
            for (i, device_name) in
                enumerate(OutputParameters["OutputData"]["DynamicDevices"])
                state_symbol = Symbol(OutputParameters["OutputData"]["States"][i])
                push!(column_names, string(device_name, "_", state_symbol))
                timeseries_data = hcat(
                    timeseries_data,
                    get_state_series(results, (device_name, state_symbol))[2],
                )
            end
        end
        d = Dict(
            "fault" => fault,
            "data" => DataFrame(timeseries_data, Symbol.(column_names)),
        )

        if OutputParameters["PVS"]["SavePVSData"]
            @assert length(pvs_data) == 4
            @assert (pvs_bus_number == bus_from) || (pvs_bus_number == bus_to)
            d["pvs"] = pvs_data
        end
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
        fault = BranchImpedanceChange(
            t_fault,
            Line,
            fault_device_name,
            fault_type[2]["Multiplier"],
        )

    elseif fault_type[1] == "ControlReferenceChange"
        g = get_component(DynamicInjection, system, fault_device_name)
        starting_Pref = get_P_ref(g)
        fault = ControlReferenceChange(
            t_fault,
            g,
            :P_ref,
            starting_Pref * fault_type[2]["RefValue"],
        )

    elseif fault_type[1] == "GeneratorTrip"
        g = get_component(DynamicInjection, system, fault_device_name)
        fault = GeneratorTrip(t_fault, g)

    elseif fault_type[1] == "LoadChange"
        l = get_component(ElectricLoad, system, fault_device_name)
        starting_Pref = get_P_ref(l)
        fault = PowerSimulationsDynamics.LoadChange(
            t_fault,
            l,
            :P_ref,
            starting_Pref * fault_type[2]["RefValue"],
        )

    elseif fault_type[1] == "LoadTrip"
        l = get_component(ElectricLoad, system, fault_device_name)
        fault = LoadTrip(t_fault, l)

    else
        error("This type of fault is not supported")
    end
    return fault
end
