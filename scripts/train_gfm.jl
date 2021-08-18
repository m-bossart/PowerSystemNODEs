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
using Flux.Losses: mae, mse
using ForwardDiff
using Statistics
using GalacticOptim
include("../models/DynamicComponents.jl")
include("../models/InverterModels.jl")
include("../models/utils.jl")
include("../models/parameter_utils.jl")
include("../models/init_functions.jl")
configure_logging(console_level = Logging.Info)

########################PARAMETERS##############################################
label = "" #base label for training figures
const train_split = 0.99
solver = Rodas4() #Rodas5()
abstol = 1e-8
reltol = 1e-5
tfault =  0.01
tspan = (0.0, 1.0)
steps = 100
tsteps =  10 .^ (range(log10(tfault), log10(tspan[2]),length= steps))
#tsteps = tspan[1]:step:tspan[2]


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
        abstol = abstol,
        reltol = reltol,
        reset_simulation=true, saveat=tsteps);

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
Ir_scale = maximum(ir_truth[2]) - minimum(ir_truth[2])
Ii_scale = maximum(ii_truth[2]) - minimum(ii_truth[2])
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

######INITIALIZE THE GFM SURROGATE AND BUILD THE TRAINING PROBLEM###############
f = get_init_gfm(p_ode, x₀[5], x₀[19])
res = nlsolve(f, x₀)
@assert converged(res)
dx = similar(x₀)
gfm(dx,res.zero,p_ode,0.0)
@assert all(isapprox.(dx, 0.0; atol=1e-8))

M = MassMatrix(19, 0)
gfm_func = ODEFunction(gfm, mass_matrix = M)
gfm_prob = ODEProblem(gfm_func,res.zero,tspan,p_ode)

sol = solve(gfm_prob, solver, abstol=abstol, reltol=reltol, saveat=tsteps)
scatter!(p3, sol, vars = [5], markersize= 2,  label = "real current gfm surrogate")
scatter!(p4, sol, vars = [19], markersize= 2, label = "imag current gfm surrogate")
display(plot(p1,p2,p3,p4, layout = (2,2), xlim = (0.0,0.1)))
##

global u₀ = res.zero
global refs

function predict_gfm(θ)
    p = vcat(θ, refs,  p_fixed)
    _prob = remake(gfm_prob, p=p, u0=u₀)
    solve(_prob,solver, abstol=abstol, reltol=reltol, saveat=tsteps)
end

function loss_gfm(θ)
    solution = predict_gfm(θ)
    sol = Array(solution)
    if solution.retcode == :Success
        loss = mae(sol[5,:], ir_truth[2] ) / Ir_scale + mae( sol[19,:], ii_truth[2]) / Ii_scale
    else
        loss = Inf
    end
    loss, sol
end

list_plots = []
list_losses = Float64[]
list_gradnorm = Float64[]
cb_gfm = function(θ,l, sol)
    #DISPLAY LOSS AND PLOT
    grad_norm = Statistics.norm(ForwardDiff.gradient(x -> first(loss_gfm(x)),θ),2)
    push!(list_gradnorm, grad_norm)
    push!(list_losses, l)
    display(l)
    cb_gfm_plot(sol)

    #UPDATE REFERENCES AND INITIAL CONDITIONS
    x₀, refs_int = initialize_sys!(sys_init, "gen1", θ)
    global refs = refs_int
    p = vcat(θ, refs, p_fixed)
    f = get_init_gfm(p, x₀[5], x₀[19])
    res = nlsolve(f, x₀)
    @assert converged(res)
    global u₀ = res.zero

    return false
end

#TRAIN THE PLAIN GFM
res_gfm = @time DiffEqFlux.sciml_train(loss_gfm, p_inv, ADAM(0.01), GalacticOptim.AutoZygote(), cb = cb_gfm, maxiters = 2) #TODO make maxiters = (#oftrainingfaults x iters_per_fault)

##
anim = Animation()
for plt in list_plots
    frame(anim, plt)
end
display(anim)
gif(anim, "figs/gfm_train.gif", fps = 3)
scatter(list_losses, title = "loss")
scatter(list_gradnorm, title = "norm of gradient")




## Not yet used below

#TODO - methods for evaluating the trained surrogates on the full test set.



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
