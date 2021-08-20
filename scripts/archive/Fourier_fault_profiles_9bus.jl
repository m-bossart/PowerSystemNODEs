using Pkg
Pkg.activate(".")

using Revise
using OrdinaryDiffEq
using PowerSystems
using PowerSimulationsDynamics
const PSID = PowerSimulationsDynamics
const PSY = PowerSystems
using Plots
using FFTW
using Statistics
include("../../models/DynamicComponents.jl")
include("../../models/InverterModels.jl")

#SIMULATION PARAMETERS
dtmax = 0.02
tspan = (0.0, 10.0)
step = 1e-2
tsteps = tspan[1]:step:tspan[2]
tfault = 0.1
solver =Rodas5()
#Parameters for reference change
fault_inj = 2
P_ref =  1.6
#Parameters for ybus change
fault_branch = "1"
r_factor = 0.5
x_factor = 0.5

sys = System("cases/case9.m")     #build system from case file
surrogate_device_name = "gen-1"   #The name of the device we want to surrogatize

list_disturbances = []

for g in  get_components(Generator, sys)
    add_component!(sys, inv_case78(get_name(g)), g)  #add gfm inverters
    #add_component!(sys, dyn_gen_classic(get_name(g)), g) #add sgs
end
surrogate_device = get_component(DynamicInverter, sys, "gen-1")

branches = deepcopy(collect(get_components(Branch, sys)))
for fault_branch in branches
    #Define a network switch disturbance
    branches = deepcopy(collect(get_components(Branch, sys)))
    for b in branches
        if b.name == fault_branch.name
            b.r = b.r*2.0
            b.x = b.x*2.0
        end
    end
    Ybus_fault = Ybus(branches, get_components(Bus, sys))[:, :]
    disturbance_NetworkSwitch = NetworkSwitch(
        tfault, #tfault
        Ybus_fault, #new Ybus
    )
    push!(list_disturbances, disturbance_NetworkSwitch)
    push!(list_disturbances, BranchTrip(tfault,fault_branch.name))
end

sim = Simulation!(
    MassMatrixModel,
    sys,
    pwd(),
    tspan,

)

#Define a control reference change disturbance
injs = collect(get_components(DynamicInjection, sys))
for fault_inj in injs

    disturbance_ControlReferenceChange = ControlReferenceChange(tfault, fault_inj , PowerSimulationsDynamics.P_ref_index,  get_ext(fault_inj)["control_refs"][PowerSimulationsDynamics.P_ref_index]*0.8)

    push!(list_disturbances, disturbance_ControlReferenceChange)
end

p1 =plot()
p2 =plot()
p3 =plot()
p4 =plot()


for dist in list_disturbances
    if typeof(dist) == NetworkSwitch
        color = :black
    elseif typeof(dist) == BranchTrip
        color =:red
    elseif typeof(dist) == ControlReferenceChange
        color = :green
    end
    sim = Simulation!(
        MassMatrixModel,
        sys,
        pwd(),
        tspan,
        dist,
    )

    execute!(sim, #simulation structure
            solver, #IDA() is Sundials DAE Solver for implicit form
            reset_simulation=true,dtmax=dtmax,saveat=tsteps); #

    surrogate_bus = get_component(Generator, sim.sys, surrogate_device_name).bus.number
    V = get_voltage_magnitude_series(sim,surrogate_bus)[2]
    θ = get_voltage_angle_series(sim,surrogate_bus)[2]
    t = sim.solution.t[1:end-1]

    plot!(p1, t, V, title = "Voltage magnitude time series",xlabel="time(s)",ylabel="V(pu)",linecolor=color)
    plot!(p2, t, θ, title = "Voltage angle time series",xlabel="time(s)",ylabel="θ(rad)",linecolor=color)

    N = length(V)
    fs = (N-1)/(tspan[2]-tspan[1])
    freqs = fftfreq(N, fs)
    freqs_pos = freqs[freqs .>= 0]
    F_V = fft(V)
    F_θ = fft(θ)

    F_V = F_V[freqs .>= 0]
    F_V = F_V/N
    F_V[2:end]= F_V[2:end]*2

    plot!(p3, freqs_pos[2:end], abs.(F_V[2:end]),title = "Voltage magnitude frequency spectrum", xlabel="freq(Hz)",ylabel="|Fᵥ|",linewidth=1,linecolor=color)
    F_θ = F_θ[freqs .>= 0]
    F_θ = F_θ/N
    F_θ[2:end]= F_θ[2:end]*2

    plot!(p4, freqs_pos[2:end], abs.(F_θ[2:end]),title="Voltage angle frequency spectrum",xlabel="freq(Hz)",ylabel="|Fθ|",linewidth=1,linecolor=color)

end
p = plot(p1,p2,p3,p4,layout=(2,2),legend=false,size =(3000,2000))
png(p,"figs/fault_fft_invs")
display(p)
