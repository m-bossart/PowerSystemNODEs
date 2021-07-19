using Pkg
Pkg.activate(".")
using HDF5
using Revise
using OrdinaryDiffEq
using PowerSystems
using PowerSimulationsDynamics
const PSID = PowerSimulationsDynamics
const PSY = PowerSystems
using Plots
using FFTW
using Statistics
include("../models/DynamicComponents.jl")
include("../models/InverterModels.jl")

#TODO Run a larger set at end of day.


#After each fault, use this function to reset to the default system and then make the modifications for that specific fault.
function base_system()
    sys = System("cases/case9.m")
    add_component!(sys, inv_case78(surrogate_device_name), get_component(Generator,sys,surrogate_device_name))
    return sys
end

function build_disturbances(sys)
    disturbances = []

    #BRANCH FAULTS
    branches = deepcopy(collect(get_components(Branch, sys)))
    for b in branches
        push!(disturbances, BranchTrip(tfault,get_name(b)))
    end

    #REFERENCE CHANGE FAULTS
    injs = collect(get_components(DynamicInjection, sys))
    injs = filter(x -> x !== get_component(DynamicInjection,sys,surrogate_device_name), injs) # Don't allow change in the surrogate injector (filter it!)

    for fault_inj in injs
        for Pref in Prefchange
            disturbance_ControlReferenceChange = ControlReferenceChange(tfault, fault_inj , PowerSimulationsDynamics.P_ref_index,  get_ext(fault_inj)["control_refs"][PowerSimulationsDynamics.P_ref_index] * Pref)
            push!(disturbances, disturbance_ControlReferenceChange)
        end
    end
    return disturbances
end

#SIMULATION PARAMETERS
dtmax = 1e-2
tspan = (0.0, 5.0)
step = 1e-3
tsteps = tspan[1]:step:tspan[2]
N = length(tsteps) #for dft
fs = (N-1)/(tspan[2]-tspan[1])
freqs = fftfreq(N, fs)
freqs_pos = freqs[freqs .>= 0]

tfault = 0.1
solver = Rodas5()
surrogate_device_name = "gen-2"   #The name of the device we want to surrogatize
devices = [inv_case78]#[inv_case78 dyn_gen_classic] # TODO add in the GFL and VSM
dispatches = [1.0]#[0.9,1.0]  #fraction of loading relative to nominal case read from file
Prefchange = [1.0] # [0.8, 0.9] #fraction of initial reference

global sys = base_system()
surrogate_bus =  get_component(Generator, sys, surrogate_device_name).bus
surrogate_bus_number = surrogate_bus.number
base_sys = base_system()

p1 =plot()
p2 =plot()
fault_sources = PeriodicVariableSource[]
counter = 0
count_stable = 0
count_unstable = 0
for a in devices #Generator at bus 15
    for b in devices #Generator at bus 16
        gens = collect(get_components(Generator, sys))
        gens = filter(x -> x !== get_component(Generator,sys,surrogate_device_name), gens) #fILTER OUT THE SURROGATE DEVICE
        dyn_models = [a, b]
        for (i,gen) in enumerate(gens)
            gen_name = get_name(gen)
            add_component!(sys, dyn_models[i](gen_name),gen)
        end
        for e in dispatches #Dispatch level of the system
            for load in get_components(PowerLoad,sys)
                orig_load = get_component(PowerLoad,base_sys,get_name(load))
                set_active_power!(load, get_active_power(orig_load)*e)
            end

            #Need to solve power flow before modifying reference values.
            sim = Simulation!(
                MassMatrixModel,
                sys,
                pwd(),
                tspan,
            )

            disturbances = []
            disturbances = build_disturbances(sys)

            for (n,f) in enumerate(disturbances)
                sim = Simulation!(
                    MassMatrixModel,
                    sys,
                    pwd(),
                    tspan,
                    f,
                )

                execute!(sim, #simulation structure
                        solver, #IDA() is Sundials DAE Solver for implicit form
                        reset_simulation=true,dtmax=dtmax,saveat=tsteps); #

                V = get_voltage_magnitude_series(sim,surrogate_bus_number)[2]
                θ = get_voltage_angle_series(sim,surrogate_bus_number)[2]
                t = sim.solution.t[1:end-1]
                if (sim.solution.retcode == :Success)
                    plot!(p1, t, V, title = "Voltage magnitude time series",xlabel="time(s)",ylabel="V(pu)",color=:black, linewidth=1,size =(3000,2000))
                    plot!(p2, t, θ, title = "Voltage angle time series",xlabel="time(s)",ylabel="θ(rad)",color=:black, linewidth=1, size =(3000,2000))

                    F_V = fft(V)
                    F_V = F_V[freqs .>= 0]
                    F_V = F_V/N
                    F_V[2:end]= F_V[2:end]*2
                    internal_voltage_coefficients = [(-imag(f), real(f)) for f in F_V[2:end]]

                    F_θ = fft(θ)
                    F_θ = F_θ[freqs .>= 0]
                    F_θ = F_θ/N
                    F_θ[2:end]= F_θ[2:end]*2
                    internal_angle_coefficients = [(-imag(f), real(f)) for f in F_V[2:end]]

                    fault_source = PeriodicVariableSource(
                        name = "InfBus",
                        available = false,
                        bus = surrogate_bus,
                        R_th = 0.0,
                        X_th = 0.0,
                        internal_voltage_bias = abs(F_V[1]),
                        internal_voltage_frequencies = freqs_pos[2:end],
                        internal_voltage_coefficients = internal_voltage_coefficients, # [( -imag.(F_V[2:end]) , real.(F_V[2:end]) )] ,
                        internal_angle_bias =  abs(F_θ[1]),
                        internal_angle_frequencies =  freqs_pos[2:end],
                        internal_angle_coefficients =internal_angle_coefficients ,
                    )

                    push!(fault_sources, fault_source)


                    count_stable += 1
                else
                    count_unstable += 1
                end
                #TODO: do the fft, write the source coefficients to PeriodicVariableSource, push to a list of sources.. To check we can write a script to loop through the script and re-generate the time domain data .
                counter  += 1;  print("Plotting simulation number ", counter, "\n")

            end

        end
        sys = base_system()

    end
end

#SUMMARIZE THE RUNS...
p = plot(p1,p2, layout= (2,1))
display(p)
png(p,"figs/fault_fft_invs")
print("Considered...\n", length(devices), " device model(s) at each of the non-surrogate generator bus\n")
print(length(dispatches), " system loading levels\n")
print(length(disturbances), " line or reference change faults\n\n\n")
print(count_stable, " stable runs\n")
print(count_unstable, " unstable runs\n")


#Check that you are
p3 = plot()
for fault_source in fault_sources
    V_reconstruct = zeros(length(tsteps))
    V_reconstruct = V_reconstruct .+ get_internal_voltage_bias(fault_source)
    freqs = get_internal_voltage_frequencies(fault_source)
    coeffs = get_internal_voltage_coefficients(fault_source)

    for (i,f) in enumerate(freqs)
        global V_reconstruct
        V_reconstruct += coeffs[i][1]* sin.(f .*2 .* pi .* tsteps)
        V_reconstruct += coeffs[i][2]* cos.(f .*2 .* pi .* tsteps)
    end
    plot!(p3,tsteps, V_reconstruct, title="reconstructed voltage mag")
end
display(p3)


sys2 = System("cases/omib.m")
add_source_to_ref(sys2)
source = collect(get_components(Source,sys2))[1]
set_dynamic_injector!(source, fault_sources[1])
add_component!(sys2, fault_sources[1], source)



## STORING THE DATASET TO DISK?
