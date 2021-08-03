#Script for exploring techniques to initialize surrogate models.

using Pkg
Pkg.activate(".")
using Revise
using NLsolve
using PowerSystems
const PSY = PowerSystems
using PowerSimulationsDynamics
using OrdinaryDiffEq
using DiffEqFlux
const PSID = PowerSimulationsDynamics

include("../models/DynamicComponents.jl")
include("../models/init_functions.jl")
# Place holders for time_dependent functions from Periodic Source
Vr(t) = 1.0000000201046026
Vi(t) = 1.0880616805155173e-9


#SIMULATION PARAMETERS
dtmax = 0.02
tspan = (0.0, 1.0)
step = Float32(0.01)
tsteps = tspan[1]:step:tspan[2]

solver =Rodas5()

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

p_inv = vcat(ω_lp,kp_pll,ki_pll,Ta,kd,kω,kq,ωf,kpv,kiv,kffv,rv,lv,kpc,kic,kffi,
              ωad,kad,lf,rf,cf,lg,rg,Vref,ωref,Pref,Qref)


#CHECK THE INITILAIZATION OF THE STANDARD GFM FROM PSID
x₀_dict = get_initial_conditions(sim)[surrogate_device_name]
x₀ =[value for (key,value) in x₀_dict]
dx = similar(x₀)
gfm(dx,x₀,p_inv,0)
@assert all(isapprox.(dx, 0.0; atol=1e-6))
f = get_init_gfm(p_inv, x₀[5], x₀[19])
res = nlsolve(f, x₀)
gfm(dx,res.zero,p_inv,0)
@assert all(isapprox.(dx, 0.0; atol=1e-8))

##
#INITIALIZE THE SURROGATE WITH NN THAT DEPENDS ON STATES
nn_states = FastChain(FastDense(length(x₀), 3, tanh),
                       FastDense(3, 2))
p_nn_states = initial_params(nn_states)
n_weights_nn_states = length(p_nn_states)
p_all_states = vcat(p_nn_states, p_inv)
x₀_nn_states = vcat(x₀, 0.0,0.0,nn_states(x₀,p_nn_states)[1],nn_states(x₀, p_nn_states)[2], x₀[5], x₀[19])
dx = similar(x₀_nn_states)
g = get_init_gfm_nn_states(p_all_states, x₀[5], x₀[19])
res_nn_states= nlsolve(g,x₀_nn_states)
@assert converged(res_nn_states)
gfm_nn_states(dx,res_nn_states.zero,p_all_states,0)
@assert all(isapprox.(dx, 0.0; atol=1e-8))

##
#INITIALIZE THE SURROGATE WITH NN THAT DEPENDS ON VOLTAGE
nn_voltage = FastChain(FastDense(2, 3, tanh),
                       FastDense(3, 2))
p_nn_voltage = initial_params(nn_voltage)
n_weights_nn_voltage = length(p_nn_voltage)
p_all_voltage = vcat(p_nn_voltage, p_inv)
x₀_nn_voltage = vcat(x₀, 0.0,0.0,nn_voltage([Vr(0),Vi(0)],p_nn_voltage)[1],nn_voltage([Vr(0),Vi(0)], p_nn_voltage)[2], x₀[5], x₀[19])
dx = similar(x₀_nn_voltage)
g = get_init_gfm_nn_voltage(p_all_voltage, x₀[5], x₀[19])
res_nn_voltage= nlsolve(g,x₀_nn_voltage)
@assert converged(res_nn_voltage)
gfm_nn_voltage(dx,res_nn_voltage.zero,p_all_voltage,0)
@assert all(isapprox.(dx, 0.0; atol=1e-8))
