#activate the environment
cd(@__DIR__)
using Pkg
Pkg.activate(".")

using Revise
using OrdinaryDiffEq
using PowerSystems
using PowerSimulationsDynamics
const PSID = PowerSimulationsDynamics
const PSY = PowerSystems
using Plots
include("DynamicComponents.jl")
include("gfm_inverter.jl")

#SIMULATION PARAMETERS
dtmax = 0.02
tspan = (0.0, .5)
step = 1e-5
tsteps = tspan[1]:step:tspan[2]
tfault = 0.1

#Parameters for reference change
fault_inj = 2
P_ref =  1.6
#Parameters for ybus change
fault_branch = "1"
r_factor = 0.5
x_factor = 0.5

sys = System("cases/case9.m")     #build system from case file
surrogate_device_name = "gen-1"   #The name of the device we want to surrogatize

for g in  get_components(Generator, sys)
    case_inv = inv_case78(get_name(g))
    add_component!(sys, case_inv, g)
end
surrogate_device = get_component(DynamicInverter, sys, "gen-1")

#Define a network switch disturbance
branches = deepcopy(collect(get_components(Branch, sys)))
for b in branches
    b.r = b.r*2.0
    b.x = b.x*2.0
end
Ybus_fault = Ybus(branches, get_components(Bus, sys))[:, :]
disturbance_NetworkSwitch = NetworkSwitch(
    tfault, #tfault
    Ybus_fault, #new Ybus
)

#Define a control reference change disturbance
injs = collect(get_components(DynamicInjection, sys))
injs_names =  [get_name(i) for i in injs]
injs = injs[sortperm(injs_names)]
disturbance_ControlReferenceChange = ControlReferenceChange(tfault, injs[1] , PowerSimulationsDynamics.P_ref_index, 0.8)

#GET PARAMETERS TO USE IN UODE SURROGATE
#pll
ω_lp = get_ω_lp(surrogate_device.freq_estimator)
kp_pll = get_kp_pll(surrogate_device.freq_estimator)
ki_pll = get_ki_pll(surrogate_device.freq_estimator)
#outer control
Ta = PSY.get_Ta(surrogate_device.outer_control.active_power)
kd = PSY.get_kd(surrogate_device.outer_control.active_power)
kω = PSY.get_kω(surrogate_device.outer_control.active_power)
kq = PSY.get_kq(surrogate_device.outer_control.reactive_power)
ωf = PSY.get_ωf(surrogate_device.outer_control.reactive_power)
#inner control
kpv   = get_kpv(surrogate_device.inner_control)
kiv   = get_kiv(surrogate_device.inner_control)
kffv   = get_kffv(surrogate_device.inner_control)
rv   = get_rv(surrogate_device.inner_control)
lv   = get_lv(surrogate_device.inner_control)
kpc   = get_kpc(surrogate_device.inner_control)
kic   = get_kic(surrogate_device.inner_control)
kffi   = get_kffi(surrogate_device.inner_control)
ωad   =get_ωad(surrogate_device.inner_control)
kad   = get_kad(surrogate_device.inner_control)
#lcl
lf = get_lf(surrogate_device.filter)
rf = get_rf(surrogate_device.filter)
cf = get_cf(surrogate_device.filter)
lg = get_lg(surrogate_device.filter)
rg = get_rg(surrogate_device.filter)

sim = Simulation!(
    MassMatrixModel,
    sys,
    pwd(),
    tspan,
    disturbance_NetworkSwitch,
)

#References
Vref = PSY.get_ext(surrogate_device)["control_refs"][1]
ωref = PSY.get_ext(surrogate_device)["control_refs"][2]
Pref = PSY.get_ext(surrogate_device)["control_refs"][3]
Qref = PSY.get_ext(surrogate_device)["control_refs"][4]

p_inv = vcat(ω_lp,kp_pll,ki_pll,Ta,kd,kω,kq,ωf,kpv,kiv,kffv,rv,lv,kpc,kic,kffi,
              ωad,kad,lf,rf,cf,lg,rg,Vref,ωref,Pref,Qref)


surrogate_bus = get_component(Generator, sim.sys, surrogate_device_name).bus.number

x₀_dict = get_initial_conditions(sim)[surrogate_device_name]
x₀ = [value for (key,value) in x₀_dict ]

vr_index  = get(PSID.get_lookup(sim.simulation_inputs), surrogate_bus, 0)
vi_index = vr_index + PSID.get_bus_count(sim.simulation_inputs)
ir_index = PSID.get_global_index(sim.simulation_inputs)[surrogate_device_name][:ir_filter]
ii_index = PSID.get_global_index(sim.simulation_inputs)[surrogate_device_name][:ii_filter]

execute!(sim, #simulation structure
        Rodas5(), #IDA() is Sundials DAE Solver for implicit form
        reset_simulation=true,dtmax=dtmax,saveat=tsteps); #

Vr(t) = sim.solution(t)[vr_index]
Vi(t) = sim.solution(t)[vi_index]
Ir(t) = sim.solution(t)[ir_index]
Ii(t) = sim.solution(t)[ii_index]

push!(x₀, Vr(0.0), Vi(0.0))

p1 = plot(get_state_series(sim,(surrogate_device_name,:ir_filter)), label = "real current PSID")
p2 = plot(get_state_series(sim,(surrogate_device_name,:ii_filter)), label = "imaginary current PSID")
p3 = plot(get_voltage_magnitude_series(sim,surrogate_bus), label = "voltage mag PSID")
p4 = plot(get_voltage_angle_series(sim,surrogate_bus), label = "voltage ang PSID")
#uses Vr(t) and Vi(t)

#Build Neural network
#dim_input =  length(x₀) #same dimenson as the states of inverter
#dim_hidden = 10
#dim_output = 1 # number of NN terms incorporated into the equations.
#nn = FastChain(FastDense(dim_input, dim_hidden, tanh),
#                       FastDense(dim_hidden, dim_hidden, tanh),
#                       FastDense(dim_hidden, dim_output))

#Build parameter vector
#p_nn= initial_params(nn)
#n_weights = length(p_nn)
#p_all = vcat(p_nn,p_inv)
#one of the key problems...finding initial conditions of a surrogate
#Train on no disturbance, and see if you can come up with an initialized


M = zeros(length(x₀) ,length(x₀) )
for i = 1:length(x₀)-2
  M[i,i] = 1.0f0
end

uodefunc = ODEFunction(uode_surrogate, mass_matrix = M)
prob_uode = ODEProblem(uodefunc,x₀,tspan,p_inv)
sol = solve(prob_uode,Rodas5(), dtmax=dtmax,saveat=tsteps)      #  Use same SOVLER!


plot!(p1, sol, vars = (0, 5), label="real current DiffEq") #5,9 - cnv real current
plot!(p2, sol, vars = (0, 19), label="imaginary current DiffEq") #19, 18 - cnv imag current
Vmag = sqrt.(Vr.(sol.t).^2 + Vi.(sol.t).^2)
Vang = tan.(Vi.(sol.t) ./ Vr.(sol.t))
plot!(p3, sol.t, Vmag, label = "voltage mag DiffEq ")
plot!(p4, sol.t, Vang, label = "voltage ang DiffEq")
display(plot(p1,p2,p3,p4,layout = (2,2),legend=true,))

#plot(get_voltage_magnitude_series(sim,surrogate_bus)[2]-Vmag)
#maximum(get_voltage_magnitude_series(sim,surrogate_bus)[2]-Vmag)
#plot(sol.t,get_state_series(sim,(surrogate_device_name,:ir_filter))[2] .- sol[5,:])
