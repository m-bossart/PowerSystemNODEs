
using Pkg
Pkg.activate(".")
using Revise
using Distributions
using OrdinaryDiffEq
using DifferentialEquations
using PowerSystems
using Logging
using PowerSimulationsDynamics
const PSID = PowerSimulationsDynamics
const PSY = PowerSystems
using Plots
using FFTW
using Statistics
using NLsolve
using DiffEqFlux
using ForwardDiff
include("../models/DynamicComponents.jl")
include("../models/InverterModels.jl")
include("../models/utils.jl")
include("../models/parameter_utils.jl")
include("../models/init_functions.jl")
configure_logging(console_level = Logging.Error)

########################PARAMETERS##############################################
const train_split = 0.99
solver = Rodas5()
dtmax = 0.001
tspan = (0.0, 1.0)
step = 0.05
tsteps = tspan[1]:step:tspan[2]
################BUILD THE TRAINING SYSTEMS FOR GENERATING TRUTH DAT#############
sys_faults = System("systems/fault_library.json")
sys_full = System("systems/base_system.json")
sys_train, sys_test = build_train_test(sys_faults, sys_full, 2, train_split, add_pvs = true)
@info "training set size:", length(collect(get_components(Source,sys_train)))
@info "test set size:", length(collect(get_components(Source,sys_test)))
to_json(sys_train,"systems/sys_train.json", force = true)
to_json(sys_test,"systems/sys_test.json", force = true)

##############################GENERATE TRUE SOLUTION ###########################
available_source = activate_next_source!(sys_train)
set_bus_from_source(available_source) #Bus voltage is used in power flow, not source voltage. Need to set bus voltage from soure internal voltage

sim = Simulation!(
    MassMatrixModel,
    sys_train,
    pwd(),
    tspan,
)

@info "train system power flow", solve_powerflow(sys_train)["flow_results"]
@info "train system power flow", solve_powerflow(sys_train)["bus_results"]

p_inv = [500.0, 0.084, 4.69, 2.0, 400.0, 20.0,0.2,1000.0,0.59,  736.0, 0.0, 0.0, 0.2,  1.27, 14.3, 0.0, 50.0,  0.2,  0.08, 0.003, 0.074, 0.2,0.01]

sys_init = build_sys_init(sys_train)
transformer = collect(get_components(Transformer2W,sys_init))[1]
pvs = collect(get_components(PeriodicVariableSource, sys_init))[1]

x₀, refs = initialize_sys!(sys_init, "gen1", p_inv)

execute!(sim,
        solver,
        reset_simulation=true,dtmax=dtmax,saveat=tsteps);

active_source = collect(get_components(Source, sys_train,  x -> PSY.get_available(x)))[1]
#V, θ = Source_to_function_of_time(active_source)
V, θ = Source_to_function_of_time(get_dynamic_injector(active_source))
p_ode = vcat(p_inv, refs, get_x(transformer) + get_X_th(pvs), get_r(transformer)+ get_R_th(pvs))

Vmag_bus = get_voltage_magnitude_series(sim, 2)
θ_bus = get_voltage_angle_series(sim, 2)
#Vmag_internal = get_state_series(sim, ("source1",:Vt))
#θ_internal = get_state_series(sim, ("source1",:θt))
#p1 , p2 = plot_pvs(tsteps, get_dynamic_injector(active_source))
p1 = plot()
p2= plot()
plot!(p1, Vmag_bus, label="bus voltage")
#plot!(p1,Vmag_internal,  label="internal voltage")
display(p1)
p2 = plot!(p2, θ_bus, label="bus angle")
#plot!(p2, θ_internal, label="internal angle")
display(plot(p1,p2, layout = (2,1)))


ir_truth, ii_truth = get_total_current_series(sim) #TODO Better to measure current at the PVS (implement method after PVS is complete)
p3 = plot(ir_truth, label = "real current true")
p4 = plot(ii_truth, label = "imag current true")
display(plot(p3,p4, layout = (2,1)))
# INITIALIZE THE PLAIN GFM
f = get_init_gfm(p_ode, x₀[5], x₀[19])
res = nlsolve(f, x₀)
@assert converged(res)
dx = similar(x₀)
gfm(dx,res.zero,p_ode,0.0)
@assert all(isapprox.(dx, 0.0; atol=1e-8))

M = MassMatrix(19, 0)
gfm_func = ODEFunction(gfm, mass_matrix = M)
gfm_prob = ODEProblem(gfm_func,res.zero,tspan,p_ode)

sol = solve(gfm_prob, Rodas5(), dtmax=dtmax, saveat=tsteps)
plot!(p3, sol, vars = [5],  label = "real current gfm surrogate")
plot!(p4, sol, vars = [19],  label = "imag current gfm surrogate")
display(plot(p3,p4, layout = (2,1)))
## INITIALIZE THE SURROGATE WITH NN THAT DEPENDS ON STATES
nn_states = FastChain(FastDense(length(x₀), 3, tanh),
                       FastDense(3, 2))
p_nn_states = initial_params(nn_states)
n_weights_nn_states = length(p_nn_states)
p_all_states = vcat(p_nn_states, p_ode)
x₀_nn_states = vcat(x₀, 0.0,0.0,nn_states(x₀,p_nn_states)[1],nn_states(x₀, p_nn_states)[2], x₀[5], x₀[19])

g = get_init_gfm_nn_states(p_all_states, x₀[5], x₀[19])
res_nn_states= nlsolve(g,x₀_nn_states)
@assert converged(res_nn_states)
dx = similar(x₀_nn_states)
gfm_nn_states(dx,res_nn_states.zero,p_all_states,0.0)
@assert all(isapprox.(dx, 0.0; atol=1e-8))

M = MassMatrix(23, 2)
gfm_nn_states_func = ODEFunction(gfm_nn_states, mass_matrix = M)
gfm_nn_states_prob = ODEProblem(gfm_nn_states_func,res_nn_states.zero,tspan,p_all_states)

sol = solve(gfm_nn_states_prob, Rodas5(), dtmax=dtmax, saveat=tsteps)
plot!(p3, sol, vars = [24],  label = "real current gfm+nn surrogate")
plot!(p4, sol, vars = [25],  label = "imag current gfm+nn surrogate")
display(plot(p3,p4, layout = (2,1)))

## INITIALIZE THE SURROGATE WITH NN THAT DEPENDS ON VOLTAGE - DON"T CONSIDER YET
#nn_voltage = FastChain(FastDense(2, 3, tanh),
#                       FastDense(3, 2))
#p_nn_voltage = initial_params(nn_voltage)
#n_weights_nn_voltage = length(p_nn_voltage)
#p_all_voltage = vcat(p_nn_voltage, p_ode)
#x₀_nn_voltage = vcat(x₀, 0.0,0.0,nn_voltage([V(0.0), θ(0.0)],p_nn_voltage)[1], nn_voltage([V(0.0), θ(0.0)],p_nn_voltage)[2], x₀[5], x₀[19])
#h = get_init_gfm_nn_voltage(p_all_voltage, x₀[5], x₀[19])
#res_nn_voltage= nlsolve(h,x₀_nn_voltage)
#@assert converged(res_nn_voltage)
#dx = similar(x₀_nn_voltage)
#gfm_nn_voltage(dx,res_nn_voltage.zero,p_all_voltage,0.0)
#@assert all(isapprox.(dx, 0.0; atol=1e-8))
#M = MassMatrix(23, 2)
#gfm_nn_voltage_func = ODEFunction(gfm_nn_voltage, mass_matrix = M)
#gfm_nn_voltage_prob = ODEProblem(gfm_nn_voltage_func,res_nn_voltage.zero,tspan,p_all_voltage)
#sol = solve(gfm_nn_voltage_prob, Rodas5(), dtmax=dtmax, saveat=tsteps)
#plot(sol,vars = [24,25],  title="GFM+NN(voltage) surrogate")
#TODO Improve the plotting of power system results (PowerGraphics?)
#TODO Make sure the surrogate matches the PSID system when we have identical models...
##

global u₀ = res.zero

function predict_gfm(θ)                                                         #TODO Constrain parameters? Or just use smaller step in optimizer?
    p = vcat(θ, refs, get_x(transformer), get_r(transformer))
    #display("predict")
    #display(u₀)
    _prob = remake(gfm_prob, p=p, u0=u₀)
    Array(solve(_prob,Rodas5(), dtmax=dtmax, saveat=tsteps))
end

function loss_gfm(θ)
    sol = predict_gfm(θ)
    loss = sum(abs2,  ir_truth[2] - sol[5,:]) + sum(abs2, ii_truth[2] - sol[19,:])
end
l = loss_gfm(p_inv)

cb_gfm = function(θ,l)
    push!(list_losses,l)
    display(l)
    x₀, refs = initialize_sys!(sys_init, "gen1", θ)
    p = vcat(θ, refs, get_x(transformer), get_r(transformer))
    f = get_init_gfm(p, x₀[5], x₀[19])
    res = nlsolve(f, x₀)
    @assert converged(res)
    global u₀ = res.zero
    #display("callback")
    #display(u₀)
    p1 = plot(predict_gfm(θ)[5,:], label = "real current prediction")
    plot!(p1, ir_truth[2], label = "real current true")
    p2 = plot(predict_gfm(θ)[19,:], label = "imag current prediction")
    plot!(p2, ii_truth[2], label = "imag current true")
    plt = plot(p1,p2,layout=(2,1))
    push!(list_plots, plt)
    display(plt)
    return false
end

#TRAIN THE PLAIN GFM
training_res = @time DiffEqFlux.sciml_train(loss_gfm, p_inv, ADAM(0.005), cb = cb_gfm, maxiters = 10) #TODO make maxiters = (#oftrainingfaults x iters_per_fault)


##
global u₀_nn_states = res_nn_states.zero
function predict_gfm_nn_states(θ)
    p = vcat(θ, p_ode)
    _prob = remake(gfm_nn_states_prob, p=p,  u0=u₀_nn_states)
    Array(solve(_prob,Rodas5(), dtmax=dtmax, saveat=tsteps, sensealg = QuadratureAdjoint(autojacvec=ReverseDiffVJP(true))))
end

function loss_gfm_nn_states(θ)
    sol = predict_gfm_nn_states(θ)
    loss = sum(abs2,  ir_truth[2] - sol[5,:]) + sum(abs2, ii_truth[2] - sol[19,:])
end


list_plots = []
list_losses = Float64[]
cb_nn_states = function(θ,l)
    push!(list_losses,l)
    display(l)
    x₀, refs = initialize_sys!(sys_init, "gen1", θ)
    p = vcat(θ, p_ode)
    x₀_nn_states = vcat(x₀, 0.0,0.0,nn_states(x₀,θ)[1],nn_states(x₀, θ)[2], x₀[5], x₀[19])
    f = get_init_gfm_nn_states(p, x₀[5], x₀[19])
    res = nlsolve(f, x₀_nn_states)
    @assert converged(res)
    global u₀_nn_states = res.zero
    display("callback")
    display(u₀)
    p1 = plot(predict_gfm_nn_states(θ)[24,:], label = "real current prediction")
    plot!(p1, ir_truth[2], label = "real current true")
    p2 = plot(predict_gfm_nn_states(θ)[25,:], label = "imag current prediction")
    plot!(p2, ii_truth[2], label = "imag current true")
    plt = plot(p1,p2,layout=(2,1))
    push!(list_plots, plt)
    display(plt)
    return false
end

trainig_res3 = @time DiffEqFlux.sciml_train(loss_gfm_nn_states, p_nn_states, ADAM(0.005), cb = cb_nn_states, maxiters = 2)


## Not yet used below
























cb_iters = function (θ,l)
    push!(list_losses,l)    #record loss
    if (loss < ϵ) # or if you've reach some number of epochs?
        activate_next_source!(sys)
        sim = Simulation!(MassMatrixModel, sys, pwd(),tspan)
        exectute!(sim)
        #TODO: Modify the ground truth data here
        active_pvs = get_components(PeriodicVariableSource, sys,  x -> PSY.get_available(x))
        V, θ = PVS_to_function_of_time(active_pvs)
    end
    return false
end
cb_threshold = function (θ,l)
    push!(list_losses,l)    #record loss

    if (loss < ϵ) # or if you've reach some number of epochs?
        activate_next_source!(sys)
        sim = Simulation!(MassMatrixModel, sys, pwd(),tspan)
        execture!(sim)
        #TODO: Modify the ground truth data here
        active_pvs = get_components(PeriodicVariableSource, sys,  x -> PSY.get_available(x))
        V, θ = PVS_to_function_of_time(active_pvs)
    end
    return false
end
