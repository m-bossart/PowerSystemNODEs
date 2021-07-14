using NLsolve
using PowerSystems
const PSY = PowerSystems
using PowerSimulationsDynamics
using OrdinaryDiffEq
const PSID = PowerSimulationsDynamics

include("../models/DynamicComponents.jl")
include("../models/init_functions.jl")

#SIMULATION PARAMETERS
dtmax = 0.02
tspan = (0.0, 1.0)
step = Float32(0.01)
tsteps = tspan[1]:step:tspan[2]

solver =Rodas5()#, #Rodas4P2() # Rodas5(), Rodas4P2(),

sys = System("cases/case9.m")     #build system from case file
surrogate_device_name = "gen-1"   #The name of the device we want to surrogatize (ie get initial conditions for)

for g in  get_components(Generator, sys)
    case_inv = inv_case78(get_name(g))
    add_component!(sys, case_inv, g)
end
surrogate_device = get_component(DynamicInverter, sys, "gen-1")

sim = Simulation!(
    MassMatrixModel,
    sys,
    pwd(),
    tspan,
)
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

#References (must come after building simulation)
Vref = PSY.get_ext(surrogate_device)["control_refs"][1]
ωref = PSY.get_ext(surrogate_device)["control_refs"][2]
Pref = PSY.get_ext(surrogate_device)["control_refs"][3]
Qref = PSY.get_ext(surrogate_device)["control_refs"][4]

p_inv = Float32.(vcat(ω_lp,kp_pll,ki_pll,Ta,kd,kω,kq,ωf,kpv,kiv,kffv,rv,lv,kpc,kic,kffi,
              ωad,kad,lf,rf,cf,lg,rg,Vref,ωref,Pref,Qref))


surrogate_bus = get_component(Generator, sim.sys, surrogate_device_name).bus.number
Vm₀ = get_initial_conditions(sim)["Vm"][surrogate_bus]
θ₀ =  get_initial_conditions(sim)["θ"][surrogate_bus]
Vr₀ = Vm₀ * cos(θ₀)
Vi₀ =  Vm₀ * sin(θ₀)
x₀_dict = get_initial_conditions(sim)[surrogate_device_name]
x₀ = Float32.([value for (key,value) in x₀_dict])

f = get_init_gfm(p_inv, x₀[1], x₀[10], x₀[9], x₀[19])

res = nlsolve(f, x₀)
