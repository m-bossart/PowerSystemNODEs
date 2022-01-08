# TODO: Use Dict temporarily during dev while the fields are defined
# WHich things change with fault? 
struct NODETrainInputs      #could move common data to fields outside of dict 
    tsteps::Vector{Float64}
    fault_data::Dict{String, Dict{Symbol, Any}}
end

function serialize(inputs::NODETrainInputs, file_path::String)
    open(file_path, "w") do io
        JSON3.write(io, inputs)
    end
    return
end

mutable struct NODETrainDataParams
    solver::String
    solver_tols::Tuple{Float64, Float64}
    tspan::Tuple{Float64, Float64}
    steps::Int64
    tsteps_spacing::String
    ode_model::String #"vsm" or "none"
    base_path::String
    output_data_path::String
end

StructTypes.StructType(::Type{NODETrainDataParams}) = StructTypes.Struct()
StructTypes.StructType(::Type{NODETrainInputs}) = StructTypes.Struct()
function NODETrainDataParams(;
    solver = "Rodas4",
    solver_tols = (1e-6, 1e-9),
    tspan = (0.0, 3.0),
    steps = 300,
    tsteps_spacing = "linear",
    ode_model = "vsm",
    base_path = pwd(),
    output_data_path = joinpath(base_path, "input_data"),
)
    NODETrainDataParams(
        solver,
        solver_tols,
        tspan,
        steps,
        tsteps_spacing,
        ode_model,
        base_path,
        output_data_path,
    )
end

#TODO - break up this function for ease of understanding/debugging 
function generate_train_data(sys_train, NODETrainDataParams, SURROGATE_BUS)
    tspan = NODETrainDataParams.tspan
    steps = NODETrainDataParams.steps
    if NODETrainDataParams.tsteps_spacing == "linear"
        tsteps = tspan[1]:((tspan[2] - tspan[1]) / steps):tspan[2]
    end
    solver = instantiate_solver(NODETrainDataParams)
    abstol = NODETrainDataParams.solver_tols[1]
    reltol = NODETrainDataParams.solver_tols[2]
    fault_data = Dict{String, Dict{Symbol, Any}}()

    for pvs in collect(get_components(PeriodicVariableSource, sys_train))
        available_source = activate_next_source!(sys_train)
        set_bus_from_source(available_source) #TODO - Bus voltage is used in power flow, not source voltage. Need to set bus voltage from soure internal voltage

        sim_full = Simulation!(
            MassMatrixModel,
            sys_train,
            pwd(),
            tspan,
            console_level = PSID_CONSOLE_LEVEL,
            file_level = PSID_FILE_LEVEL,
        )

        execute!(
            sim_full,
            solver,
            abstol = abstol,
            reltol = reltol,
            reset_simulation = false,
            saveat = tsteps,
        )
        active_source =
            collect(get_components(Source, sys_train, x -> PSY.get_available(x)))[1]
        ode_data = get_total_current_series(sim_full)

        transformer = collect(get_components(Transformer2W, sys_train))[1]
        p_network = [get_x(transformer) + get_X_th(pvs), get_r(transformer) + get_R_th(pvs)]

        #flow_results = solve_powerflow(sys_train)["flow_results"]
        bus_results = solve_powerflow(sys_train)["bus_results"]
        @info "full system", bus_results
        surrogate_bus_result = bus_results[in([SURROGATE_BUS]).(bus_results.bus_number), :]
        #=         if flow_results[1,:bus_from] == SURROGATE_BUS
                    P_pf = flow_results[1,:P_from_to]      
                    Q_pf = flow_results[1,:Q_from_to]
                elseif flow_results[1,:bus_to] == SURROGATE_BUS
                    P_pf = flow_results[1,:P_to_from]     
                    Q_pf = flow_results[1,:Q_to_from]
                else 
                    @error "Surrogate bus not found in the system"
                end  =#
        P_pf = surrogate_bus_result[1, :P_gen] / get_base_power(sys_train)   #TODO- divide by base power  
        Q_pf = surrogate_bus_result[1, :Q_gen] / get_base_power(sys_train)
        V_pf = surrogate_bus_result[1, :Vm]
        θ_pf = surrogate_bus_result[1, :θ]

        @info collect(get_components(Bus, sys_train))
        @warn "P*", P_pf, "Q*", Q_pf, "V*", V_pf, "θ*", θ_pf

        fault_data[get_name(pvs)] = Dict(
            :p_network => p_network,
            :p_pf => [P_pf, Q_pf, V_pf, θ_pf],  #Expand to include all PF data [P,Q,V,θ] (p_pf)
            :ir_ground_truth => ode_data[1, :],
            :ii_ground_truth => ode_data[2, :],
        )
        if NODETrainDataParams.ode_model == "vsm"
            #################### BUILD INITIALIZATION SYSTEM ###############################
            sys_init, p_inv = build_sys_init(sys_train) #returns p_inv, the set of average parameters 
            x₀, refs, Vr0, Vi0 = initialize_sys!(sys_init, "gen1")
            Vm, Vθ = Source_to_function_of_time(get_dynamic_injector(active_source))
            p_ode = vcat(p_inv, refs)
            sim_simp = Simulation!(
                MassMatrixModel,
                sys_init,
                pwd(),
                tspan,
                console_level = PSID_CONSOLE_LEVEL,
                file_level = PSID_FILE_LEVEL,
            )
            @info "initialize system power flow", solve_powerflow(sys_init)["flow_results"]
            @info "initialize system power flow", solve_powerflow(sys_init)["bus_results"]
            @debug show_states_initial_value(sim_simp)
            @time execute!(
                sim_simp,
                solver,
                abstol = abstol,
                reltol = reltol,
                initializealg = NoInit(),
                reset_simulation = false,
                saveat = tsteps,
            )

            avgmodel_data = get_total_current_series(sim_simp)
            @warn get_name(pvs)

            fault_data[get_name(pvs)][:p_ode] = p_ode
            fault_data[get_name(pvs)][:x₀] = x₀
            fault_data[get_name(pvs)][:ir_node_off] = avgmodel_data[1, :]
            fault_data[get_name(pvs)][:ii_node_off] = avgmodel_data[2, :]
        end
    end
    return NODETrainInputs(tsteps, fault_data)
end
