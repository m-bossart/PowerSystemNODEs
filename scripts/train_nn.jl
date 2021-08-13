#TODO Improve the plotting of power system results (PowerGraphics?)
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
configure_logging(console_level = Logging.Info)

########################PARAMETERS##############################################
const train_split = 0.99
solver = Rodas4() #Rodas5()
dtmax = 0.001
tspan = (0.0, 1.0)
step = 0.02
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
execute!(sim,
        solver,
        reset_simulation=true, dtmax=dtmax, saveat=tsteps);

active_source = collect(get_components(Source, sys_train,  x -> PSY.get_available(x)))[1]
Vmag_bus = get_voltage_magnitude_series(sim, 2)
θ_bus = get_voltage_angle_series(sim, 2)
Vmag_internal = get_state_series(sim, ("source1",:Vt))
θ_internal = get_state_series(sim, ("source1",:θt))
p1 , p2 = plot_pvs(tsteps, get_dynamic_injector(active_source))
plot!(p1, Vmag_bus, label="bus voltage")
plot!(p1,Vmag_internal,  label="internal voltage")
p2 = plot!(p2, θ_bus, label="bus angle")
plot!(p2, θ_internal, label="internal angle")
ir_truth, ii_truth = get_total_current_series(sim) #TODO Better to measure current at the PVS (implement method after PVS is complete)
p3 = plot(ir_truth, label = "real current true")
p4 = plot(ii_truth, label = "imag current true")

#####################BUILD INITIALIZATION SYSTEM################################
p_inv = [500.0, 0.084, 4.69, 2.0, 400.0, 20.0,0.2,1000.0,0.59,  736.0, 0.0, 0.0, 0.2,  1.27, 14.3, 0.0, 50.0,  0.2,  0.08, 0.003, 0.074, 0.2,0.01]
sys_init = build_sys_init(sys_train)
transformer = collect(get_components(Transformer2W,sys_init))[1]
pvs = collect(get_components(PeriodicVariableSource, sys_init))[1]
p_fixed =  [get_x(transformer) + get_X_th(pvs), get_r(transformer)+ get_R_th(pvs)]
x₀, refs = initialize_sys!(sys_init, "gen1", p_inv)
V, θ = Source_to_function_of_time(get_dynamic_injector(active_source))
p_ode = vcat(p_inv, refs, p_fixed)

######INITIALIZE THE GFM+NN SURROGATE AND BUILD THE TRAINING PROBLEM############
nn = FastChain(FastDense(2, 3, tanh),
                       FastDense(3, 2))
p_nn = initial_params(nn)
n_weights_nn = length(p_nn)
p_all = vcat(p_nn, p_inv, refs, p_fixed, n_weights_nn)
x₀_nn = vcat(x₀, 0.0,0.0,nn([V(0.0), θ(0.0)],p_nn)[1], nn([V(0.0), θ(0.0)],p_nn)[2], x₀[5], x₀[19])
h = get_init_gfm_nn(p_all, x₀[5], x₀[19])
res_nn= nlsolve(h, x₀_nn)
@assert converged(res_nn)
dx = similar(x₀_nn)
gfm_nn(dx,res_nn.zero,p_all,0.0)
@assert all(isapprox.(dx, 0.0; atol=1e-8))
M = MassMatrix(23, 2)
gfm_nn_func = ODEFunction(gfm_nn, mass_matrix = M)
gfm_nn_prob = ODEProblem(gfm_nn_func,res_nn.zero,tspan,p_all)
sol = solve(gfm_nn_prob, solver, dtmax=dtmax, saveat=tsteps)
sol = solve(gfm_nn_prob, solver, dtmax=dtmax, saveat=tsteps)
plot!(p3, sol, vars = [24],  label = "real current gfm+nn surrogate")
plot!(p4, sol, vars = [25],  label = "imag current gfm+nn surrogate")
display(plot(p1,p2,p3,p4, layout = (2,2)))

##

global u₀ = res_nn.zero
global refs

function predict_gfm_nn(θ)
    p = vcat(θ, p_inv, refs, p_fixed, n_weights_nn)
    _prob = remake(gfm_nn_prob, p=p,  u0=u₀)
    Array(solve(_prob, solver, dtmax=dtmax, saveat=tsteps))
end

function loss_gfm_nn(θ)
    sol = predict_gfm_nn(θ)
    loss = sum(abs2,  ir_truth[2] - sol[5,:]) + sum(abs2, ii_truth[2] - sol[19,:])
    loss, sol
end

list_plots = []
list_losses = Float64[]
cb_gfm_nn = function(θθ,l, sol)
    #DISPLAY LOSS AND PLOT
    push!(list_losses,l)
    display(l)
    cb_gfm_nn_plot(sol)

    #UPDATE REFERENCES AND INITIAL CONDITIONS
    x₀, refs_int = initialize_sys!(sys_init, "gen1", p_inv) #TODO don't need this each time
    global refs = refs_int
    p = vcat(θθ, p_inv, refs, p_fixed, n_weights_nn)
    x₀_nn = vcat(x₀, 0.0, 0.0, nn([V(0.0), θ(0.0)],θθ)[1], nn([V(0.0), θ(0.0)],θθ)[2], x₀[5], x₀[19])
    f = get_init_gfm_nn(p, x₀[5], x₀[19])
    res = nlsolve(f, x₀_nn)
    @assert converged(res)
    global u₀ = res.zero

    return false
end

trainig_res3 = @time DiffEqFlux.sciml_train(loss_gfm_nn, p_nn, ADAM(0.01), cb = cb_gfm_nn, maxiters = 20)

##
anim = Animation()
for plt in list_plots
    frame(anim, plt)
end
display(anim)
gif(anim, "figs/nm_train.gif", fps = 3)
