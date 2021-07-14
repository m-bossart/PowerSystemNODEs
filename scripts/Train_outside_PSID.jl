#activate the environment
using Pkg
Pkg.activate(".")

using Revise
using DiffEqFlux
using OrdinaryDiffEq
using PowerSystems
using PowerSimulationsDynamics
const PSID = PowerSimulationsDynamics
const PSY = PowerSystems
using Plots
include("../models/DynamicComponents.jl")
include("../models/InverterModels.jl")

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


##
M = Float32.(zeros(length(x₀) ,length(x₀) ))
for i = 1:length(x₀)   #-2     Include if using the IB version (last two equations are algebraic)
  M[i,i] = 1.0f0
end

f_train = 5
A_train = 0.01
Vr(t) =  Vr₀
Vr(t) = Vr₀ + A_train * sin.(2*pi*f_train*t)
Vi(t) =  Vi₀

func_gfm = ODEFunction(gfm, mass_matrix = M)
prob_gfm = ODEProblem(func_gfm,x₀,Float32.(tspan),p_inv)
sol_gfm = solve(prob_gfm,solver, dtmax=dtmax,saveat=tsteps)      #  Use same SOVLER!


p1 = plot( sol_gfm, vars = (0, 5), label="real current DiffEq") #5,9 - cnv real current
p2 = plot( sol_gfm, vars = (0, 19), label="imaginary current DiffEq") #19, 18 - cnv imag current
Vmag = sqrt.(Vr.(sol_gfm.t).^2 + Vi.(sol_gfm.t).^2)
Vang = tan.(Vi.(sol_gfm.t) ./ Vr.(sol_gfm.t))
p3 = plot(sol_gfm.t, Vmag, label = "voltage mag DiffEq ")
p4 = plot(sol_gfm.t, Vang, label = "voltage ang DiffEq")
p5 = plot(sol_gfm)
p6 = plot(sol_gfm)


#BUILD NN AND UDE PROBLEM

dim_input =  length(x₀)
dim_hidden = 3
dim_output = 2
nn = FastChain(FastDense(dim_input, dim_hidden, tanh),
                       FastDense(dim_hidden, dim_output))
p_nn= initial_params(nn)
n_weights = length(p_nn)
p_inv_mod = p_inv .*  (1 .+ ((rand(length(p_inv)).-0.5)./5))  #parameters are either +/- 10%
#p_all = vcat(p_nn,p_inv_mod)
p_all = vcat(p_nn, p_inv)
func_gfm_nn_dist = ODEFunction(gfm_nn_dist, mass_matrix = M)
prob_gfm_nn_dist = ODEProblem(func_gfm_nn_dist,x₀,Float32.(tspan),p_all)
sol_gfm_nn_dist = solve(prob_gfm_nn_dist,solver, dtmax=dtmax,saveat=tsteps)


#COMPARE FULL MODEL AND UODE BEFORE TRAINING
plot!(p1, sol_gfm_nn_dist, vars = (0, 5), label="real current uode") #5,9 - cnv real current
plot!(p2, sol_gfm_nn_dist, vars = (0, 19), label="imaginary current uode") #19, 18 - cnv imag current
Vmag = sqrt.(Vr.(sol_gfm_nn_dist.t).^2 + Vi.(sol_gfm_nn_dist.t).^2)
Vang = tan.(Vi.(sol_gfm_nn_dist.t) ./ Vr.(sol_gfm_nn_dist.t))
plot!(p3, sol_gfm_nn_dist.t, Vmag, label = "voltage mag uode ")
plot!(p4, sol_gfm_nn_dist.t, Vang, label = "voltage ang uode")
plot!(p5, sol_gfm_nn_dist[1:10,:])
plot!(p6, sol_gfm_nn_dist[11:19,:])
display(plot(p1,p2,p3,p4,p5,p6,layout = (3,2),legend=true))


#TRAINING
function predict_solution(θ)
    p = vcat(θ,p_inv)   #p_inv is fixed, only the NN parameters are updated
    #p = θ
    Array(solve(prob_gfm_nn_dist,Rodas5(),p=p, saveat = tsteps, sensealg=InterpolatingAdjoint())) #autojacvec=ReverseDiffVJP(true)
end

function loss(θ)
 pred = predict_solution(θ)
 loss = sum(abs2, ( pred[5,:] .- sol_gfm[5,:] ))     #real current part
 loss += sum(abs2, ( pred[19,:] .- sol_gfm[19,:] ))  #imaginary current part
end

list_loss = []
list_grad = []
list_plots = []

iter = 1
n_plot = 1  #plot every step
cb = function (θ,l)
    global list_loss, list_plots, iter
    push!(list_loss, l)
    print(l,"\n")
    if iter % n_plot == 0
        pred = predict_solution(θ)
        p1 = scatter(pred[5,:], label="prediction")
        p2 = scatter(pred[19,:],label="prediction")

        plot!(p1,  sol_gfm[5,:], label="full model", title = string("real current, ",string(iter)))
        plot!(p2,  sol_gfm[19,:], label="full model", title = string("imaginary current, ",string(iter)))
        plt = plot(p1,p2)
        png(plt, "figs/loss")
        push!(list_plots, plt)
        display(plt)
    end
    iter += 1
    return false
end


@time loss(p_nn)
##

res = @time DiffEqFlux.sciml_train(loss,  p_nn,ADAM(0.1), cb = cb, maxiters= 100)

anim = Animation()
for plt in list_plots
    frame(anim, plt)
end
display(anim)
gif(anim,fps = 1)
