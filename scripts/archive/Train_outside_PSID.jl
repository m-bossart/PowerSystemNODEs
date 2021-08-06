#NOTE: p_inv (inverter parameters), p_ode (all ode parameters), p_nn (nn parameters)
#TODO everything should be Float32 when you have NNs?
using Pkg
Pkg.activate(".")
using Revise
using DiffEqFlux
using OrdinaryDiffEq
using PowerSystems
using PowerSimulationsDynamics
const PSID = PowerSimulationsDynamics
const PSY = PowerSystems
using NLsolve
using Plots
include("../../models/DynamicComponents.jl")
include("../../models/InverterModels.jl")
include("../../models/utils.jl")
include("../../models/parameter_utils.jl")
include("../../models/init_functions.jl")

#SIMULATION PARAMETERS
dtmax = 0.002f0
tspan = (0.0f0, 1.0f0)
step = 0.01f0
tsteps = tspan[1]:step:tspan[2]

solver =Rodas5() #, #Rodas4P2() # Rodas5(), Rodas4P2(),


sys = System("cases/case9.m")
surrogate_device_name = "gen-1"

surrogate_bus = get_component(Generator, sys, surrogate_device_name).bus.number

for g in  get_components(Generator, sys)
    case_inv = inv_case78(get_name(g))
    add_component!(sys, case_inv, g)
end
surrogate_device = get_component(DynamicInverter, sys, surrogate_device_name)
p_inv = get_parameters(surrogate_device)
x₀ , refs = initialize_sys!(sys,surrogate_device_name,p_inv)
#TODO cleaner way to get Vm₀ and θ₀
 sim = Simulation!(
        MassMatrixModel,
        sys,
        pwd(),
        (0.0, 1.0),
    )
Vm₀ = get_initial_conditions(sim)["Vm"][surrogate_bus]
θ₀ =  get_initial_conditions(sim)["θ"][surrogate_bus]

V(t) = Vm₀ + 0.01* sin(20*t)
θ(t) = θ₀  + 0.01* sin(10*t)


##
M = MassMatrix(length(x₀), 0)


p_ode = vcat(p_inv,refs, 0.0f0, 0.0f0)

f = get_init_gfm(p_ode, x₀[5], x₀[19])

res = nlsolve(f, x₀)
@assert converged(res)

func_gfm = ODEFunction(gfm, mass_matrix = M)
prob_gfm = ODEProblem(func_gfm,Float32.(res.zero),tspan,p_ode)
sol_gfm = solve(prob_gfm,solver, dtmax=dtmax,saveat=tsteps)


p1 = plot( sol_gfm, vars = (0, 5), label="real current DiffEq") #5,9 - cnv real current
p2 = plot( sol_gfm, vars = (0, 19), label="imaginary current DiffEq") #19, 18 - cnv imag current
Vmag = V.(sol_gfm.t)
Vang = θ.(sol_gfm.t)
p3 = plot(sol_gfm.t, Vmag, label = "voltage mag DiffEq ")
p4 = plot(sol_gfm.t, Vang, label = "voltage ang DiffEq")
p5 = plot(sol_gfm)
p6 = plot(sol_gfm)


#BUILD NN AND UDE PROBLEM

dim_input =  length(x₀)
dim_hidden = 3
dim_output = 2
nn_states = FastChain(FastDense(dim_input, dim_hidden, tanh),
                       FastDense(dim_hidden, dim_output))
p_nn_states= initial_params(nn_states)
n_weights_nn_states = length(p_nn_states)
#p_inv_mod = p_inv .*  (1 .+ ((rand(length(p_inv)).-0.5)./5))  #parameters are either +/- 10%

p_all = vcat(p_nn_states, p_ode)
f = get_init_gfm_nn_states(p_all, x₀[5], x₀[19])
x₀_nn_states = vcat(x₀, 0.0,0.0,nn_states(x₀,p_nn_states)[1],nn_states(x₀, p_nn_states)[2], x₀[5], x₀[19])
res = nlsolve(f, x₀_nn_states)
@assert converged(res)

M2 = MassMatrix(23, 2)
#TODO Update gfm_nn_dist to the newer version, once this is up and running, contact JDL
func_gfm_nn_states = ODEFunction(gfm_nn_states, mass_matrix = M2)
prob_gfm_nn_states = ODEProblem(func_gfm_nn_states, Float32.(x₀_nn_states), tspan, p_all)
sol_gfm_nn_states = solve(prob_gfm_nn_states, solver, dtmax=dtmax,saveat=tsteps)


#COMPARE FULL MODEL AND UODE BEFORE TRAINING
plot!(p1, sol_gfm_nn_states, vars = (0, 24), label="real current uode") #5,9 - cnv real current
plot!(p2, sol_gfm_nn_states, vars = (0, 25), label="imaginary current uode") #19, 18 - cnv imag current
Vmag = V.(sol_gfm_nn_states.t)
Vang = θ.(sol_gfm_nn_states.t)
plot!(p3, sol_gfm_nn_states.t, Vmag, label = "voltage mag uode ")
plot!(p4, sol_gfm_nn_states.t, Vang, label = "voltage ang uode")
plot!(p5, sol_gfm_nn_states[1:10,:])
plot!(p6, sol_gfm_nn_states[11:19,:])
display(plot(p1,p2,p3,p4,p5,p6,layout = (3,2),legend=true))

##
#TRAINING
sensealg = InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true))
#sensealg = QuadratureAdjoint()
#sensealg = ReverseDiffAdjoint()
p_ode = Float32.(p_ode)
function predict_solution(θ)
    p = vcat(θ,p_ode)   #p_ode is fixed, only the NN parameters are updated
    print(typeof(p))
    Array(solve(prob_gfm_nn_states ,Rodas5(), p=p, saveat = tsteps, sensealg=sensealg)) #autojacvec=ReverseDiffVJP(true)
end

function loss(θ)
 pred = predict_solution(θ)
 loss = sum(abs2, ( pred[24,:] .- sol_gfm[5,:] ))     #real current part
 loss += sum(abs2, ( pred[25,:] .- sol_gfm[19,:] ))  #imaginary current part
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


@time loss(p_nn_states)
##

res = @time DiffEqFlux.sciml_train(loss,  p_nn_states, ADAM(0.1), cb = cb, maxiters= 4)

anim = Animation()
for plt in list_plots
    frame(anim, plt)
end
display(anim)
gif(anim,fps = 1)
