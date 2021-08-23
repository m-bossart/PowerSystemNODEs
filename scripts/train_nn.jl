
include("../models/DynamicComponents.jl")
include("../models/InverterModels.jl")
include("../models/utils.jl")
include("../models/parameter_utils.jl")
include("../models/init_functions.jl")
configure_logging(console_level = Logging.Error)

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

######INITIALIZE THE GFM+NN SURROGATE AND BUILD THE TRAINING PROBLEM############
nn = FastChain(FastDense(2, nn_width, nn_activation),
                FastDense(nn_width, nn_width, nn_activation),
                FastDense(nn_width, 2))

p_nn = initial_params(nn)
n_weights_nn = length(p_nn)
p_all = vcat(p_nn, p_inv, refs, p_fixed, nn([V(0.0), θ(0.0)],p_nn)[1], nn([V(0.0), θ(0.0)],p_nn)[2])
x₀_nn = vcat(x₀, 0.0, 0.0, x₀[5], x₀[19])
h = get_init_gfm_nn(p_all, x₀[5], x₀[19])
res_nn= nlsolve(h, x₀_nn)
@assert converged(res_nn)
dx = similar(x₀_nn)
gfm_nn(dx,res_nn.zero,p_all,0.0)
@assert all(isapprox.(dx, 0.0; atol=1e-8))
M = MassMatrix(21, 2)
gfm_nn_func = ODEFunction(gfm_nn, mass_matrix = M)
gfm_nn_prob = ODEProblem(gfm_nn_func,res_nn.zero,tspan,p_all)

sol = solve(gfm_nn_prob, solver,  abstol=abstol, reltol=reltol,  saveat=tsteps)
scatter!(p3, sol, vars = [22], markersize=2, label = "real current gfm+nn surrogate")
scatter!(p4, sol, vars = [23], markersize=2, label = "imag current gfm+nn surrogate")
display(plot(p1,p2,p3,p4, layout = (2,2)))

##

global u₀ = res_nn.zero
global refs

function predict_gfm_nn(θθ)
    p = vcat(θθ, p_inv, refs, p_fixed, nn([V(0.0), θ(0.0)],θθ)[1],  nn([V(0.0), θ(0.0)],θθ)[2] )
    _prob = remake(gfm_nn_prob, p=p,  u0=u₀)
    #solve(_prob, solver,  abstol=abstol, reltol=reltol, saveat=tsteps  ) #
    #solve(_prob, solver,  abstol=abstol, reltol=reltol, saveat=tsteps, sensealg = QuadratureAdjoint(autojacvec=ReverseDiffVJP(true)))
    #solve(_prob, solver,  abstol=abstol, reltol=reltol, saveat=tsteps, sensealg = QuadratureAdjoint(autojacvec=ZygoteVJP()))
    solve(_prob, solver,  abstol=abstol, reltol=reltol, saveat=tsteps, sensealg = ForwardDiffSensitivity() )
end

function loss_gfm_nn(θ)
    solution = predict_gfm_nn(θ)
    sol = Array(solution)
    #if solution.retcode == :Success
    if size(sol)[2] == size(tsteps)[1]
        loss = mae(sol[22,:], ir_truth[2] ) / Ir_scale + mae( sol[23,:], ii_truth[2]) / Ii_scale
    else
        loss = Inf
    end
    loss, sol
end

list_plots = []
list_losses = Float64[]
list_gradnorm = Float64[]
cb_gfm_nn = function(θθ,l, sol) #TODO Deal with  θθ nonsense and make better naming
    #DISPLAY LOSS AND PLOT
    grad_norm = Statistics.norm(ForwardDiff.gradient(x -> first(loss_gfm_nn(x)),θθ),2)
    push!(list_gradnorm, grad_norm)
    push!(list_losses,l)
    display(l)
    cb_gfm_nn_plot(sol)

    #UPDATE REFERENCES AND INITIAL CONDITIONS
    x₀, refs_int = initialize_sys!(sys_init, "gen1", p_inv) #TODO don't need this each time
    global refs = refs_int
    p = vcat(θθ, p_inv, refs, p_fixed, nn([V(0.0), θ(0.0)],θθ)[1], nn([V(0.0), θ(0.0)],θθ)[2] )
    x₀_nn = vcat(x₀, 0.0, 0.0,  x₀[5], x₀[19])
    f = get_init_gfm_nn(p, x₀[5], x₀[19])
    res = nlsolve(f, x₀_nn)
    @assert converged(res)
    global u₀ = res.zero

    return false
end

#TRAIN THE PLAIN GFM
res_gfm = @time DiffEqFlux.sciml_train(loss_gfm_nn, p_nn, optimizer, GalacticOptim.AutoForwardDiff(), cb = cb_gfm_nn, maxiters = 1000) #TODO make maxiters = (#oftrainingfaults x iters_per_fault)

anim = Animation()
for plt in list_plots
    frame(anim, plt)
end

gif(anim, string("figs/", label ,"_train.gif"), fps = 20)
png(plot(list_losses), string("figs/",label,"_loss.png"))
png(plot(list_gradnorm), string("figs/", label,"_gradnorm.png"))

#ForwardDiff.gradient(x -> first(loss_gfm_nn(x)), p_nn)
#t_F = @elapsed ForwardDiff.gradient(x -> first(loss_gfm_nn(x)), p_nn)
#@info "ForwardDiff time", t_F

# loss_gfm_nn(p_nn)
# t_solve = @elapsed loss_gfm_nn(p_nn)

# Zygote.gradient(x -> first(loss_gfm_nn(x)), p_nn)
# t_Z = @elapsed Zygote.gradient(x -> first(loss_gfm_nn(x)), p_nn)

#aa = @time ForwardDiff.gradient(x -> first(loss_gfm_nn(x)), p_nn)
#res_nn = @time DiffEqFlux.sciml_train(loss_gfm_nn, p_nn, ADAM(0.01), GalacticOptim.AutoForwardDiff(), cb = cb_gfm_nn, maxiters = 2)
#res_nn = @time DiffEqFlux.sciml_train(loss_gfm_nn, p_nn, ADAM(0.01), GalacticOptim.AutoZygote(), cb = cb_gfm_nn, maxiters = 2)


# loss_gfm_nn(p_nn)
# t_solve = @elapsed loss_gfm_nn(p_nn)

# Zygote.gradient(x -> first(loss_gfm_nn(x)), p_nn)
# t_Z = @elapsed Zygote.gradient(x -> first(loss_gfm_nn(x)), p_nn)


# @info "Solve time Quadrature ", t_solve
# @info "Zygote gradient time Quadrature", t_Z


# ##

# scatter(list_losses, title = "loss")
# scatter(list_gradnorm, title = "norm of gradient")
