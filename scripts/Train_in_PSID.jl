using Pkg
Pkg.activate(".")
using Revise
using Distributions
using OrdinaryDiffEq
using DifferentialEquations
using PowerSystems
using PowerSimulationsDynamics
const PSID = PowerSimulationsDynamics
const PSY = PowerSystems
using Plots
using FFTW
using Statistics
using NLsolve
using DiffEqFlux
include("../models/DynamicComponents.jl")
include("../models/InverterModels.jl")
include("../models/StaticComponents.jl")
include("../models/utils.jl")
include("../models/init_functions.jl")
########################PARAMETERS##############################################
const train_split = 0.99    #proportion of faults for training (rest for test)
n_devices = 10             #number of devices in the surrogate
param_range = (0.9,1.1)
total_rating = 150.0  #MVA rating at the surrogate bus.

solver = Rodas5()
dtmax = 0.02
tspan = (0.0, 1.0)
step = 1e-2
tsteps = tspan[1]:step:tspan[2]

################BUILD THE TRAINING SYSTEMS FOR GENERATING TRUTH DAT#############
sys_faults = System("systems/fault_library.json")
sys_full = System("systems/base_system.json")
sys_train, sys_test = build_train_test(sys_faults, sys_full, 2, train_split, add_pvs = false) #Don't add PVS because can't build a sim with it yet.
@info "training set size:", length(collect(get_components(Source,sys_train)))
@info "test set size:", length(collect(get_components(Source,sys_test)))

to_json(sys_train,"systems/sys_train.json", force = true )
to_json(sys_test,"systems/sys_test.json", force = true)

##############################TRAINING##########################################
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

for g in get_components(DynamicInverter,sys_train)
    @info "real current", get_initial_conditions(sim)[get_name(g)][:ir_filter]
    @info "imag current", get_initial_conditions(sim)[get_name(g)][:ii_filter]
end
p_inv = [500.0, 0.084, 4.69, 2.0, 400.0, 20.0,0.2,1000.0,0.59,  736.0, 0.0, 0.0, 0.2,  1.27, 14.3, 0.0, 50.0,  0.2,  0.08, 0.003, 0.074, 0.2,0.01]

sys_init = build_sys_init(sys_train)              #Build the initialization system (done once)
transformer = collect(get_components(Transformer2W,sys_init))[1]

x₀, refs = initalize_sys_init!(sys_init, p_inv) #Set the current parameters, get the initial conditions and refs


@info "init system power flw", solve_powerflow(sys_init)["flow_results"]
@info "init system power flow", solve_powerflow(sys_init)["bus_results"]

execute!(sim,
        solver,
        reset_simulation=true,dtmax=dtmax,saveat=tsteps);



active_source = collect(get_components(Source, sys_train,  x -> PSY.get_available(x)))[1]
V, θ = Source_to_function_of_time(active_source)
p_all = vcat(p_inv, refs, get_x(transformer), get_r(transformer))

ir_truth, ii_truth = get_total_current_series(sim)
plot(ir_truth, title="PSID system (truth)")
plot!(ii_truth)
## INITIALIZE THE PLAIN GFM

f = get_init_gfm(p_all, x₀[5], x₀[19])
res = nlsolve(f, x₀)
@assert converged(res)
dx = similar(x₀)
gfm(dx,res.zero,p_all,0.0)
@assert all(isapprox.(dx, 0.0; atol=1e-8))

M = MassMatrix(19, 0)
gfm_func = ODEFunction(gfm, mass_matrix = M)
gfm_prob = ODEProblem(gfm_func,res.zero,tspan,p_all)

sol = solve(gfm_prob, Rodas5(), dtmax=dtmax, saveat=tsteps)
plot(sol, vars = [5,19], title="GFM surrogate")
##

## INITIALIZE THE SURROGATE WITH NN THAT DEPENDS ON STATES

nn_states = FastChain(FastDense(length(x₀), 3, tanh),
                       FastDense(3, 2))
p_nn_states = initial_params(nn_states)
n_weights_nn_states = length(p_nn_states)
p_all_states = vcat(p_nn_states, p_all)
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
plot(sol, vars = [24,25], title="GFM+NN(states) surrogate")
## INITIALIZE THE SURROGATE WITH NN THAT DEPENDS ON VOLTAGE
nn_voltage = FastChain(FastDense(2, 3, tanh),
                       FastDense(3, 2))
p_nn_voltage = initial_params(nn_voltage)
n_weights_nn_voltage = length(p_nn_voltage)
p_all_voltage = vcat(p_nn_voltage, p_all)
x₀_nn_voltage = vcat(x₀, 0.0,0.0,nn_voltage([V(0.0), θ(0.0)],p_nn_voltage)[1], nn_voltage([V(0.0), θ(0.0)],p_nn_voltage)[2], x₀[5], x₀[19])
h = get_init_gfm_nn_voltage(p_all_voltage, x₀[5], x₀[19])
res_nn_voltage= nlsolve(h,x₀_nn_voltage)
@assert converged(res_nn_voltage)
dx = similar(x₀_nn_voltage)
gfm_nn_voltage(dx,res_nn_voltage.zero,p_all_voltage,0.0)
@assert all(isapprox.(dx, 0.0; atol=1e-8))

M = MassMatrix(23, 2)
gfm_nn_voltage_func = ODEFunction(gfm_nn_voltage, mass_matrix = M)
gfm_nn_voltage_prob = ODEProblem(gfm_nn_voltage_func,res_nn_voltage.zero,tspan,p_all_voltage)

sol = solve(gfm_nn_voltage_prob, Rodas5(), dtmax=dtmax, saveat=tsteps)
plot(sol,vars = [24,25],  title="GFM+NN(voltage) surrogate")

#TODO Improve the plotting of power system results
#TODO Will need to implement a method similar to get_activepower_series() for the PVS.
    #Will require a compute_output_current method for the PVS.
    #I can implement this once Jose has implemented the PVS dynamic model.
#TODO Make sure the surrogate matches the PSID system when we have identical models...


## PSEUDO CODE BELOW
function predict_gfm(θ)
    x₀, refs = initalize_sys_init!(sys_init, θ)
    p_all = vcat(θ, refs, get_x(transformer), get_r(transformer))
    f = get_init_gfm(p_all, x₀[5], x₀[19])
    res = nlsolve(f, x₀)
    @assert converged(res)

    gfm_prob = ODEProblem(gfm_func,res.zero,tspan,p_all)
    #u₀ = initialize_surrogate()
    _prob = remake(gfm_prob, p=p_all, u0=res.zero)
    Array(solve(_prob,Rodas5(), dtmax=dtmax, saveat=tsteps))
end

function loss_gfm(θ)
    ir_pred = predict_gfm(θ)[5,:]
    ii_pred = predict_gfm(θ)[19,:]
    loss = sum(abs2,  ir_truth[2] - ir_pred)
    loss += sum(abs2, ii_truth[2] - ii_pred)
end

l = loss_gfm(p_inv)

function predict_gfm_nn_states(θ)        #Start with the non-UODE case, just optimize parameters.
    u₀ = initialize_surrogate()
    Array(solve(prob_surrogate))
end

function loss_gfm_nn_states(θ)
    pred = predict_solution(θ)
    loss = sum(abs2,  ir_truth - ir_pred)
    loss += sum(abs2, ii_truth - ii_pred)
end

function predict_gfm_nn_voltage(θ)        #Start with the non-UODE case, just optimize parameters.
    u₀ = initialize_surrogate()
    Array(solve(prob_surrogate))
end

function loss_gfm_nn_voltage(θ)
    pred = predict_solution(θ)
    loss = sum(abs2,  ir_truth - ir_pred)
    loss += sum(abs2, ii_truth - ii_pred)
end






list_losses = Float64[]

#TODO Need to figure out how to modify funcitons/data from within the callback.
cb = function(θ,l)
    push!(list_losses,l)
    return false
end
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

res = @time DiffEqFlux.sciml_train(loss_gfm, p_inv, ADAM(0.1), cb = cb, maxiters = 5) #TODO make maxiters = (#oftrainingfaults x iters_per_fault)
res = @time DiffEqFlux.sciml_train(loss,  p_nn, ADAM(0.1), cb = cb_threshold, maxiters = 10 )
