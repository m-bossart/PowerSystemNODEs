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
include("../models/InverterModels.jl")
include("../models/StaticComponents.jl")
include("../models/utils.jl")
include("../models/init_functions.jl")
#TODO Run a larger set at end of day

function add_devices_to_surrogatize!(sys::System, n_devices::Integer, surrogate_bus_number::Integer, inf_bus_number:: Integer)
    param_range = (0.9,1.1)
    surrogate_bus = collect(get_components(Bus,sys,x->get_number(x)==surrogate_bus_number))[1]
    inf_bus = collect(get_components(Bus,sys,x->get_number(x)==inf_bus_number))[1]

    surrogate_area = Area(;name = "surrogate")
    add_component!(sys,surrogate_area)
    set_area!(surrogate_bus, surrogate_area)
    set_area!(inf_bus, surrogate_area)

    gens = collect(get_components(ThermalStandard, sys, x->get_number(get_bus(x)) == surrogate_bus_number))

    !(length(gens) == 1) && @error "number of devices at surrogate bus not equal to one"
    gen = gens[1]
    total_rating = get_rating(gen) #doesn't impact dynamics
    total_base_power = get_base_power(gen)
    total_active_power = get_active_power(gen)
    remove_component!(sys,gen)
    for i in 1:n_devices
        g = ThermalStandard(
           name = string(i),
           available = true,
           status = true,
           bus = surrogate_bus,
           active_power = total_active_power, #Only divide base power by n_devices
           reactive_power = 0.0,
           rating =  total_rating/n_devices,
           active_power_limits=(min=0.0, max=3.0),
           reactive_power_limits= (-3.0,3.0),
           ramp_limits=nothing,
           operation_cost=ThreePartCost(nothing),
           base_power =  total_base_power/n_devices,
           )
       add_component!(sys, g)
       inv_typ = inv_case78(get_name(g))
       randomize_inv_parameters!(inv_typ, param_range)
       add_component!(sys, inv_typ, g)
   end
end

function build_disturbances(sys) #TODO make this more flexible, add options for which faults to include
    disturbances = []
    #BRANCH FAULTS
    lines = deepcopy(collect(get_components(Line, sys, x-> get_name(x) == "BUS 4       -BUS 5       -i_7"))) #TODO - change back to
    for l in lines
        push!(disturbances, BranchTrip(tfault,get_name(l)))
    end
    #REFERENCE CHANGE FAULTS
    injs = collect(get_components(DynamicInjection, sys,  x -> !(get_name(x) in get_name.(surrogate_gens))))
    for fault_inj in injs
        for Pref in Prefchange
            disturbance_ControlReferenceChange = ControlReferenceChange(tfault, fault_inj , PowerSimulationsDynamics.P_ref_index,  get_ext(fault_inj)["control_refs"][PowerSimulationsDynamics.P_ref_index] * Pref)
            #push!(disturbances, disturbance_ControlReferenceChange)
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
base_system_path = "systems\\base_system.json"

#BUG possible power flow issues using 14-bus system due to Fixed Admittance
#sys = System("cases/IEEE 14 bus_modified_33.raw")
sys = System("cases/IEEE 14 bus_modified_33_RemoveFixedAdmittance.raw")
#You want a single line connecting the source_bus and surrogate_bus becaues this will be the structure in the 2bus train system.
source_bus = 2      #Bus number for the bus that will become the IB
surrogate_bus = 16  #Bus number where the devices to be surrogatized are attached (cannot be reference bus in full system)

devices = [inv_case78]#[inv_case78 dyn_gen_classic] # TODO add in the GFL and VSM
#dispatches = [1.0]#[0.9,1.0]  #fraction of loading relative to nominal case read from file
Prefchange = [1.0] # [0.8, 0.9] #fraction of initial reference
n_devices = 2

add_devices_to_surrogatize!(sys, n_devices, surrogate_bus, source_bus)
to_json(sys,base_system_path , force=true)


global sys = System(base_system_path)
base_sys = System(base_system_path)
surrogate_gens = collect(get_components(ThermalStandard, sys, x-> get_bus(x).number == surrogate_bus))  #TODO reformulate so that this is a list of the surrogate generators. Because we will need to filter often to not include this set.

source_surrogate_branch = find_acbranch(source_bus, surrogate_bus)

p1 =plot()
p2 =plot()

sys_faults = System(100.0)
Bus(
    number=0,
    name="init",
    bustype=nothing,
    angle=0.0,
    magnitude=0.0,
    voltage_limits=(min=0.0, max=0.0),
    base_voltage=nothing,
    area=nothing,
    load_zone=nothing,
    ext=Dict{String, Any}(),
)
add_component!(sys_faults,Bus(number=1, name="1", bustype = BusTypes.REF, angle=0.0, magnitude=0.0, voltage_limits=(min = 0.9, max = 1.1), base_voltage=345.0))
slack_bus = [b for b in get_components(Bus, sys_faults) if b.bustype == BusTypes.REF][1]

counter = 0
count_stable = 0
count_unstable = 0
for a in devices
    for b in devices
        for c in devices
            for d in devices
                gens = get_components(ThermalStandard, sys, x -> !(get_name(x) in get_name.(surrogate_gens)))
                dyn_models = [a, b, c, d]
                for (i,gen) in enumerate(gens)
                    gen_name = get_name(gen)
                    add_component!(sys, dyn_models[i](gen_name),gen)
                end

                sim = Simulation!(MassMatrixModel, sys, pwd(), tspan) #Need to initialize before building disturbances
                disturbances = []
                disturbances = build_disturbances(sys)
                print_device_states(sim)


                for gen in get_components(ThermalStandard,sys)
                    @info get_ext(gen)
                end
                for (n,e) in enumerate(disturbances)
                    sim = Simulation!(MassMatrixModel, sys, pwd(), tspan, e)
                    to_json(sys,"systems/full_system.json", force=true)
                    @info solve_powerflow(sys)["flow_results"]
                    @info solve_powerflow(sys)["bus_results"]

                    P = get_active_power_flow(source_surrogate_branch)
                    Q = get_reactive_power_flow(source_surrogate_branch)
                    execute!(sim,
                            solver,
                            reset_simulation=true,dtmax=dtmax,saveat=tsteps);

                    V = get_voltage_magnitude_series(sim,source_bus)[2]
                    θ = get_voltage_angle_series(sim,source_bus)[2]

                    t = sim.solution.t[1:end-1]
                    if (sim.solution.retcode == :Success)
                        global count_stable += 1
                        plot!(p1, t, V, title = "Voltage magnitude time series source bus",xlabel="time(s)",ylabel="V(pu)",color=:black, linewidth=1,size =(3000,2000))
                        plot!(p2, t, θ, title = "Voltage angle time series source bus",xlabel="time(s)",ylabel="θ(rad)",color=:black, linewidth=1, size =(3000,2000))

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

                        inf_source = Source(
                            name = string(count_stable),
                            active_power = P,
                            available = false, #availability
                            reactive_power = Q,
                            bus = slack_bus, #bus
                            R_th = 0.0, #Rth
                            X_th = 5e-9,#5e-6, #Xth
                            internal_voltage = V[1],
                            internal_angle =   θ[1],
                        )
                         @info abs(F_V[1])
                         @info abs abs(F_θ[1])

                        fault_source = PeriodicVariableSource(
                            name = get_name(inf_source),
                            R_th = 0.0,
                            X_th = 0.0,
                            internal_voltage_bias = abs(F_V[1]),
                            internal_voltage_frequencies = freqs_pos[2:end],
                            internal_voltage_coefficients = internal_voltage_coefficients,
                            internal_angle_bias = abs(F_θ[1]),
                            internal_angle_frequencies =  freqs_pos[2:end],
                            internal_angle_coefficients =internal_angle_coefficients ,
                        )

                        add_component!(sys_faults, inf_source)
                        add_component!(sys_faults, fault_source, inf_source)
                    else
                        global count_unstable += 1
                    end
                    global counter  += 1;  print("Plotting simulation number ", counter, "\n")
                end

                global sys = System(base_system_path)
                @info "Building new system"
            end
        end
    end
end

#SUMMARIZE THE RUNS...
p = plot(p1,p2, layout= (2,1))
display(p)
png(p,"figs/fault_fft_invs")
print( "Considered...\n", length(devices), " device model(s) at each of the non-surrogate generator bus\n")
print(length(dispatches), " system loading levels\n")
print(length(build_disturbances(sys)), " line or reference change faults\n\n")
print("Resulting in...\n")
print(count_stable, " stable runs\n")
print(count_unstable, " unstable runs\n")


#Check that you are handling the PVS correctly by re-constructing the time domain signal
p3 = plot()
fault_sources = get_components(PeriodicVariableSource, sys_faults)
for fault_source in fault_sources
    V_reconstruct = zeros(length(tsteps))
    V_reconstruct = V_reconstruct .+ get_internal_voltage_bias(fault_source)
    freqs = get_internal_voltage_frequencies(fault_source)
    coeffs = get_internal_voltage_coefficients(fault_source)

    for (i,f) in enumerate(freqs)

        V_reconstruct += coeffs[i][1]* sin.(f .*2 .* pi .* tsteps)
        V_reconstruct += coeffs[i][2]* cos.(f .*2 .* pi .* tsteps)
    end
    plot!(p3,tsteps, V_reconstruct, title="reconstructed voltage mag")
end
display(p3)





to_json(sys_faults,"systems/fault_library.json", force=true)
