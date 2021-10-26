using Pkg
Pkg.activate(".")
using Revise
using Distributions
using OrdinaryDiffEq
using PowerSystems
using PowerSimulationsDynamics
const PSID = PowerSimulationsDynamics
const PSY = PowerSystems
using Plots
using FFTW

include("../models/DynamicComponents.jl")
include("../models/SurrogateModels.jl")
include("../models/utils.jl")
include("../models/init_functions.jl")
include("../models/parameter_utils.jl")

#SIMULATION PARAMETERS
tspan = (0.0, 2.0)  #changed from (0.0,2.0)
abstol = 1e-6
reltol = 1e-3
stepsize = 1e-3
tsteps = tspan[1]:stepsize:tspan[2]
tfault = 0.1
solver = Rodas4()

base_system_path = "systems\\base_system_3invs_vsms_20%lossP.json"
fault_system_path = "systems\\fault_library_3invs_vsms_20%lossP.json"
smooth_signal = true

#SYSTEM PARAMETERS
#BUG possible power flow issues using 14-bus system due to Fixed Admittance
#sys = System("cases/IEEE 14 bus_modified_33.raw")
sys = System("cases/IEEE 14 bus_modified_33_RemoveFixedAdmittance.raw")
source_bus = 2  #TODO Must have single line between source_bus and surrogate_bus?
surrogate_bus = 16
devices = [inv_case78]  # [inv_gfoll]# [inv_darco_droop]# [inv_case78] #  #     inv_darco_droop doesn't work, can't find the frequency?? 
Prefchange = [0.8]
n_devices = 3

add_devices_to_surrogatize!(sys, n_devices, surrogate_bus, source_bus)
to_json(sys, base_system_path, force = true)

global sys = System(base_system_path)
base_sys = System(base_system_path)

surrogate_gens =
    collect(get_components(ThermalStandard, sys, x -> get_bus(x).number == surrogate_bus))
source_surrogate_branch = find_acbranch(source_bus, surrogate_bus)

p1 = plot()
p2 = plot()

sys_faults = System(100.0)
add_component!(
    sys_faults,
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
slack_bus = [b for b in get_components(Bus, sys_faults) if b.bustype == BusTypes.REF][1]

counter = 0
count_stable = 0
count_unstable = 0
for a in devices
    for b in devices
        for c in devices
            for d in devices
                gens = get_components(
                    ThermalStandard,
                    sys,
                    x -> !(get_name(x) in get_name.(surrogate_gens)),
                )
                dyn_models = [a, b, c, d]
                for (i, gen) in enumerate(gens)
                    gen_name = get_name(gen)
                    println("GENERATOR NAMES", gen_name)
                    add_component!(sys, dyn_models[i](gen_name), gen)
                end
                #to_json(sys, "14bus_testsystem.json", force=true )
                sim = Simulation!(MassMatrixModel, sys, pwd(), tspan) #Need to initialize before building disturbances
                disturbances = []
                disturbances = build_disturbances(sys)
                println(disturbances)
                for (n, e) in enumerate(disturbances)
                    sim = Simulation!(MassMatrixModel, sys, pwd(), tspan, e)
                    to_json(sys, "systems/full_system.json", force = true)
                    @info solve_powerflow(sys)["flow_results"]
                    @info solve_powerflow(sys)["bus_results"]

                    P = get_active_power_flow(source_surrogate_branch)
                    Q = get_reactive_power_flow(source_surrogate_branch)
                    execute!(
                        sim,
                        solver,
                        abstol = abstol,
                        reltol = reltol,
                        reset_simulation = false,
                        saveat = tsteps,
                    )
                    sol = read_results(sim)
                    V_org = get_voltage_magnitude_series(sol, source_bus)[2]
                    θ_org = get_voltage_angle_series(sol, source_bus)[2]

                    tsteps_org = tsteps
                    tsteps, V = add_tanh(tsteps_org, V_org)
                    tsteps, θ = add_tanh(tsteps_org, θ_org)
                    @show length(tsteps)
                    t = sol.solution.t[1:(end - 1)]
                    if (sol.solution.retcode == :Success)
                        global count_stable += 1
                        N = length(tsteps) #for dft
                        fs = (N - 1) / (tsteps[end] - tsteps[1])
                        freqs = fftfreq(N, fs)
                        freqs_pos = freqs[freqs .>= 0] * (2 * pi)
                        plot!(
                            p1,
                            t,
                            V_org,
                            title = "Voltage magnitude time series source bus",
                            xlabel = "time(s)",
                            ylabel = "V(pu)",
                            color = :black,
                            linewidth = 1,
                        )
                        plot!(
                            p2,
                            t,
                            θ_org,
                            title = "Voltage angle time series source bus",
                            xlabel = "time(s)",
                            ylabel = "θ(rad)",
                            color = :black,
                            linewidth = 1,
                        )
                        @show fs
                        @show freqs
                        @show length(freqs_pos)
                        F_V = fft(V)
                        @show length(F_V)
                        @show length([freqs .>= 0][1])
                        @show length(V)
                        F_V = F_V[freqs .>= 0]
                        F_V = F_V / N
                        F_V[2:end] = F_V[2:end] * 2  #/ (2*pi)
                        internal_voltage_coefficients =
                            [(-imag(f), real(f)) for f in F_V[2:end]]

                        F_θ = fft(θ)
                        F_θ = F_θ[freqs .>= 0]
                        F_θ = F_θ / N
                        F_θ[2:end] = F_θ[2:end] * 2 #/ (2*pi)
                        internal_angle_coefficients =
                            [(-imag(f), real(f)) for f in F_θ[2:end]]

                        inf_source = Source(
                            name = string("source", string(count_stable)),
                            active_power = P,
                            available = false, #availability
                            reactive_power = Q,
                            bus = slack_bus, #bus
                            R_th = 0.0, #Rth
                            X_th = 5e-6, #Xth
                            internal_voltage = V[1],
                            internal_angle = θ[1],
                        )
                        println("dc voltage offset:", F_V[1])
                        println("dc angle offset:", F_θ[1])
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

                        add_component!(sys_faults, inf_source)
                        add_component!(sys_faults, fault_source, inf_source)
                    else
                        global count_unstable += 1
                    end
                    global counter += 1
                    print("Plotting simulation number ", counter, "\n")
                end

                global sys = System(base_system_path)
                @info "Building new system"
            end
        end
    end
end

#Display and Summarize Time Domain Signals
p = plot(p1, p2, layout = (2, 1))
display(p)
png(p, "figs/fault_fft_invs")
print(
    "Considered...\n",
    length(devices),
    " device model(s) at each of the non-surrogate generator bus\n",
)
print(length(build_disturbances(sys)), " line or reference change faults\n\n")
print("Resulting in...\n")
print(count_stable, " stable runs\n")
print(count_unstable, " unstable runs\n")

if false  #Check that reconstructed time domain signal matches
    #t_reconstruct= 0:stepsize:(tspan[2]-tfault)   #CHANGED
    p3 = plot()
    p4 = plot()
    fault_sources = get_components(PeriodicVariableSource, sys_faults)
    for fault_source in fault_sources
        t_reconstruct = tsteps#tspan[1]:(stepsize*10):tspan[2]
        V_reconstruct = zeros(length(t_reconstruct))
        dV_reconstruct = zeros(length(t_reconstruct))
        V_reconstruct = V_reconstruct .+ get_internal_voltage_bias(fault_source)
        retrieved_freqs = get_internal_voltage_frequencies(fault_source)
        coeffs = get_internal_voltage_coefficients(fault_source)
        for (i, ω) in enumerate(retrieved_freqs)
            V_reconstruct += coeffs[i][1] * sin.(ω .* t_reconstruct)
            V_reconstruct += coeffs[i][2] * cos.(ω .* t_reconstruct)
            dV_reconstruct +=
                ω * coeffs[i][1] * cos.(ω .* t_reconstruct) -
                ω * coeffs[i][2] * sin.(ω .* t_reconstruct)
        end
        plot!(p3, t_reconstruct, V_reconstruct, title = "reconstructed voltage mag")
        plot!(
            p4,
            t_reconstruct,
            dV_reconstruct,
            title = "voltage mag derivative",
            ylim = (-1, 1),
        )
    end
    display(plot(p3, p4, layout = (2, 1)))
end

to_json(sys_faults, fault_system_path, force = true)
