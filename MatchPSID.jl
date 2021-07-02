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

#SIMULATION PARAMETERS
dtmax = 2e-5
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


#We want to have V,I as a function of t, not a time series of discrete points because the UDE solve won't hit the same time points.

surrogate_bus = get_component(Generator, sim.sys, surrogate_device_name).bus.number

x₀_dict = get_initial_conditions(sim)[surrogate_device_name]
x₀ = [value for (key,value) in x₀_dict ]

vr_index  = get(PSID.get_lookup(sim.simulation_inputs), surrogate_bus, 0)
vi_index = vr_index + PSID.get_bus_count(sim.simulation_inputs)
ir_index = PSID.get_global_index(sim.simulation_inputs)[surrogate_device_name][:ir_filter]
ii_index = PSID.get_global_index(sim.simulation_inputs)[surrogate_device_name][:ii_filter]

execute!(sim, #simulation structure
        Rodas5(), #IDA() is Sundials DAE Solver for implicit form
        dtmax =dtmax, reset_simulation=true); #

Vr(t) = sim.solution(t)[vr_index]
Vi(t) = sim.solution(t)[vi_index]
Ir(t) = sim.solution(t)[ir_index]
Ii(t) = sim.solution(t)[ii_index]

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
function uode_surrogate(dx,x,p,t)
    #STATE INDEX AND STATES
    i__vi_filter, vi_filter = 1, x[1]
    i__γd_ic, γd_ic = 2, x[2]
    i__vq_pll, vq_pll = 3, x[3]
    i__γq_ic, γq_ic = 4, x[4]
    i__ir_filter, ir_filter = 5, x[5]
    i__ξd_ic, ξd_ic = 6, x[6]
    i__ϕd_ic, ϕd_ic = 7, x[7]
    i__ε_pll, ε_pll = 8, x[8]
    i__ir_cnv, ir_cnv = 9, x[9]
    i__vr_filter, vr_filter = 10, x[10]
    i__ω_oc, ω_oc = 11, x[11]
    i__ξq_ic, ξq_ic = 12, x[12]
    i__vd_pll, vd_pll = 13, x[13]
    i__q_oc, q_oc = 14, x[14]
    i__ϕq_ic, ϕq_ic = 15, x[15]
    i__θ_pll, θ_pll = 16, x[16]
    i__θ_oc, θ_oc = 17, x[17]
    i__ii_cnv, ii_cnv = 18, x[18]
    i__ii_filter, ii_filter = 19, x[19]

    #PARAMETERS
    ω_lp = p[1]
    kp_pll = p[2]
    ki_pll = p[3]
    Ta = p[4]
    kd = p[5]
    kω = p[6]
    kq = p[7]
    ωf = p[8]
    kpv = p[9]
    kiv = p[10]
    kffv = p[11]
    rv = p[12]
    lv = p[13]
    kpc = p[14]
    kic = p[15]
    kffi = p[16]
    ωad = p[17]
    kad = p[18]
    lf = p[19]
    rf = p[20]
    cf = p[21]
    lg = p[22]
    rg = p[23]
    Vref = p[24]    #Reference treated as parameter
    ωref = p[25]    #Reference treated as parameter
    Pref = p[26]    #Reference treated as parameter
    Qref = p[27]    #Reference treated as parameter

    ω_base = 60.0*2*pi
    ω_sys = 1.0

    #PLL
    δω_pll = kp_pll * atan(vq_pll/vd_pll) + ki_pll* ε_pll
    ω_pll = δω_pll + ω_sys
    vd_filt_pll = sin(θ_pll + pi/2)*vr_filter - cos(θ_pll + pi/2)*vi_filter
    vq_filt_pll = cos(θ_pll + pi/2)*vr_filter + sin(θ_pll + pi/2)*vi_filter

    dx[i__vd_pll]=  ω_lp*(vd_filt_pll - vd_pll)              #docs:(1a)
    dx[i__vq_pll]=  ω_lp*(vq_filt_pll - vq_pll)              #docs:(1b)
    dx[i__ε_pll] =  atan(vq_pll/vd_pll)                        #docs:(1c)
    dx[i__θ_pll] =  ω_base * δω_pll                            #docs:(1d)

    #OUTER LOOP CONTROL
    pe = vr_filter*ir_filter + vi_filter*ii_filter           #docs:(2d)
    qe = vi_filter*ir_filter - vr_filter*ii_filter           #docs:(2e)
    v_ref_olc = Vref + kq * (Qref - q_oc)


    dx[i__ω_oc]= (Pref - pe - kd*(ω_oc - ω_pll) - kω*(ω_oc - ωref)) / Ta #docs:(2a)
    dx[i__θ_oc]= ω_base*(ω_oc-ω_sys)                          #docs:(2b)
    dx[i__q_oc]= ωf * (qe - q_oc)                             #docs:(2c)

    #INNER LOOP CONTROL
    #reference transormations
    vd_filt_olc = sin(θ_oc + pi/2)*vr_filter - cos(θ_oc + pi/2)*vi_filter
    vq_filt_olc = cos(θ_oc + pi/2)*vr_filter + sin(θ_oc + pi/2)*vi_filter
    id_filt_olc = sin(θ_oc + pi/2)*ir_filter - cos(θ_oc + pi/2)*ii_filter
    iq_filt_olc = cos(θ_oc + pi/2)*ir_filter + sin(θ_oc + pi/2)*ii_filter
    id_cnv_olc =  sin(θ_oc + pi/2)*ir_cnv - cos(θ_oc + pi/2)*ii_cnv
    iq_cnv_olc =  cos(θ_oc + pi/2)*ir_cnv + sin(θ_oc + pi/2)*ii_cnv

    #Voltage control equations
    Vd_filter_ref = v_ref_olc - rv * id_filt_olc + ω_oc * lv * iq_filt_olc    #docs:(3g)
    Vq_filter_ref = -rv * iq_filt_olc - ω_oc * lv * id_filt_olc               #docs:(3h)
    dx[i__ξd_ic]=  Vd_filter_ref -  vd_filt_olc               #docs:(3a)
    dx[i__ξq_ic]=  Vq_filter_ref -  vq_filt_olc               #docs:(3b)

    #current control equations
    Id_cnv_ref =(
      kpv*(Vd_filter_ref - vd_filt_olc) + kiv*ξd_ic -                         #docs:(3i)
      cf*ω_oc*vq_filt_olc + kffi*id_filt_olc
    )
    Iq_cnv_ref =(
      kpv*(Vq_filter_ref - vq_filt_olc) + kiv*ξq_ic +                         #docs:(3j)
      cf*ω_oc*vd_filt_olc + kffi*iq_filt_olc
    )
    dx[i__γd_ic]=  Id_cnv_ref - id_cnv_olc                    #docs:(3c)
    dx[i__γq_ic]=  Iq_cnv_ref - iq_cnv_olc                    #docs:(3d)

    #active damping equations
    Vd_cnv_ref =(
      kpc*(Id_cnv_ref - id_cnv_olc) + kic*γd_ic - lf*ω_oc*iq_cnv_olc +        #docs:(3k)
      kffv*vd_filt_olc - kad*(vd_filt_olc-ϕd_ic)
    )
    Vq_cnv_ref =(
      kpc*(Iq_cnv_ref - iq_cnv_olc) + kic*γq_ic + lf*ω_oc*id_cnv_olc +        #docs:(3l)
      kffv*vq_filt_olc - kad*(vq_filt_olc-ϕq_ic)
    )
    dx[i__ϕd_ic] =  ωad * (vd_filt_olc - ϕd_ic)                #docs:(3e)
    dx[i__ϕq_ic]=  ωad * (vq_filt_olc - ϕq_ic)                #docs:(3f)

    #LCL FILTER
    #reference transformations
    Vr_cnv =   sin(θ_oc + pi/2)*Vd_cnv_ref + cos(θ_oc + pi/2)*Vq_cnv_ref
    Vi_cnv =  -cos(θ_oc + pi/2)*Vd_cnv_ref + sin(θ_oc + pi/2)*Vq_cnv_ref

    dx[i__ir_cnv]= (ω_base/lf) *                            #docs:(5a)
      (Vr_cnv - vr_filter - rf*ir_cnv + ω_sys*lf*ii_cnv)
    dx[i__ii_cnv]= (ω_base/lf) *                             #docs:(5b)
      (Vi_cnv - vi_filter - rf*ii_cnv - ω_sys*lf*ir_cnv)
    dx[i__vr_filter]= (ω_base/cf) *                       #docs:(5c)
      (ir_cnv - ir_filter + ω_sys*cf*vi_filter)
    dx[i__vi_filter]= (ω_base/cf) *                       #docs:(5d)
      (ii_cnv - ii_filter - ω_sys*cf*vr_filter)
    dx[i__ir_filter]= (ω_base/lg) *                       #docs:(5e)
      (vr_filter - Vr(t) - rg*ir_filter + ω_sys*lg*ii_filter)
    dx[i__ii_filter]= +(ω_base/lg) *                       #docs:(5f)
      (vi_filter - Vi(t) - rg*ii_filter - ω_sys*lg*ir_filter)
end

M = zeros(length(x₀) ,length(x₀) )
for i = 1:sum(length(x₀) )
  M[i,i] = 1.0f0
end

uodefunc = ODEFunction(uode_surrogate, mass_matrix = M)
prob_uode = ODEProblem(uodefunc,x₀,tspan,p_inv)
sol = solve(prob_uode,Rodas5(),p=p_inv,dtmax =dtmax)      #  Use same SOVLER!


plot!(p1, sol, vars = (0, 5), label="real current DiffEq") #5,9 - cnv real current
plot!(p2, sol, vars = (0, 19), label="imaginary current DiffEq") #19, 18 - cnv imag current
Vmag = sqrt.(Vr.(sol.t).^2 + Vi.(sol.t).^2)
Vang = tan.(Vi.(sol.t) ./ Vr.(sol.t))
plot!(p3, sol.t, Vmag, label = "voltage mag DiffEq ")
plot!(p4, sol.t, Vang, label = "voltage ang DiffEq")
savefig(plot(p1,p2,p3,p4,layout = (2,2),legend=true,xlim=(0.05,0.2)),"test.png")
length(Vmag)
length(get_voltage_magnitude_series(sim,surrogate_bus)[2])

#check all the states initial conditions
