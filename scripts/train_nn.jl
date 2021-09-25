include("../models/constants.jl")
include("../models/DynamicComponents.jl")
include("../models/InverterModels.jl")
include("../models/utils.jl")
include("../models/parameter_utils.jl")
include("../models/init_functions.jl")
configure_logging(console_level = Logging.Error)

################BUILD THE TRAINING SYSTEMS FOR GENERATING TRUTH DATA#############
sys_faults = System("systems/fault_library.json")
sys_full = System("systems/base_system.json")
sys_train, sys_test = build_train_test(sys_faults, sys_full, 2, train_split, add_pvs = true)
@info "training set size:", length(collect(get_components(Source,sys_train)))
@info "test set size:", length(collect(get_components(Source,sys_test)))
to_json(sys_train,"systems/sys_train.json", force = true)
to_json(sys_test,"systems/sys_test.json", force = true)


############################# GENERATE TRUE SOLUTION ###########################
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
println("time to solve train system for generating truth data:")
@time execute!(sim,
        solver,
        abstol = abstol,
        reltol = reltol,
        reset_simulation=true, saveat=tsteps);

active_source = collect(get_components(Source, sys_train,  x -> PSY.get_available(x)))[1]
Vmag_bus = get_voltage_magnitude_series(sim, 2)
θ_bus = get_voltage_angle_series(sim, 2)
Vmag_internal = get_state_series(sim, ("source1",:Vt))
θ_internal = get_state_series(sim, ("source1",:θt))
p1 , p2 = plot_pvs(tsteps, get_dynamic_injector(active_source), :lin)
p1_log, p2_log = plot_pvs(tsteps, get_dynamic_injector(active_source), :log)
plot!(p1, Vmag_bus, label="bus voltage")
plot!(p1,Vmag_internal,  label="internal voltage")
plot!(p1_log, Vmag_bus, label="bus voltage")
plot!(p1_log,Vmag_internal,  label="internal voltage")
plot!(p2, θ_bus, label="bus angle")
plot!(p2, θ_internal, label="internal angle")
plot!(p2_log, θ_bus, label="bus angle")
plot!(p2_log, θ_internal, label="internal angle")
#plot the other bus voltage too? to see if it is significantly different?? 
ode_data = get_total_current_series(sim) #TODO Better to measure current at the PVS (implement method after PVS is complete)
p3 = plot(tsteps, ode_data[1,:], label = "real current true")
p4 = plot(tsteps, ode_data[2,:], label = "imag current true")
p3_log = plot(tsteps,ode_data[1,:], label = "real current true", xaxis=:log)
p4_log = plot(tsteps,ode_data[2,:], label = "imag current true", xaxis=:log)

#Calculate scale factors for loss function based on true data
Ir_scale = maximum(ode_data[1,:]) - minimum(ode_data[1,:])
Ii_scale = maximum(ode_data[2,:]) - minimum(ode_data[2,:])

#################### BUILD INITIALIZATION SYSTEM ###############################
sys_init, p_inv = build_sys_init(sys_train) #returns p_inv, the set of average parameters 
transformer = collect(get_components(Transformer2W,sys_init))[1]
pvs = collect(get_components(PeriodicVariableSource, sys_init))[1]
p_fixed =  [get_x(transformer) + get_X_th(pvs), get_r(transformer)+ get_R_th(pvs)]
x₀, refs, Vr0, Vi0 = initialize_sys!(sys_init, "gen1")
Vm, Vθ = Source_to_function_of_time(get_dynamic_injector(active_source))
p_ode = vcat(p_inv, refs, p_fixed)
sim = Simulation!(
    MassMatrixModel,
    sys_init,
    pwd(),
    tspan,
)
@time execute!(sim,
        solver,
        abstol = abstol,
        reltol = reltol,
        reset_simulation=true, saveat=tsteps);
avgmodel_data = get_total_current_series(sim)

##### INITIALIZE THE GFM+NN SURROGATE AND BUILD THE TRAINING PROBLEM ###########
nn = build_nn(2, 2, nn_width, nn_hidden, nn_activation)
p_nn = initial_params(nn)
n_weights_nn = length(p_nn)
p_all = vcat(p_nn, p_inv, refs, p_fixed, Vr0, Vi0)
x₀_nn = vcat(x₀, 0.0, 0.0, x₀[5], x₀[19])

#h = get_init_gfm_nn(p_all, x₀[5], x₀[19])
#res_nn= nlsolve(h, x₀_nn)
#@assert converged(res_nn)
#dx = similar(x₀_nn)
#gfm_nn(dx,res_nn.zero,p_all,0.0)
#@assert all(isapprox.(dx, 0.0; atol=1e-8))
M = MassMatrix(21, 2)
gfm_nn_func = ODEFunction(gfm_nn, mass_matrix = M)
gfm_nn_prob = ODEProblem(gfm_nn_func, x₀_nn, tspan, p_all)

sol = solve(gfm_nn_prob, solver,  abstol=abstol, reltol=reltol,  saveat=tsteps )
scatter!(p3, sol, vars = [22], markersize=1, label = "real current gfm+nn surrogate")
scatter!(p4, sol, vars = [23], markersize=1, label = "imag current gfm+nn surrogate")
scatter!(p3_log, sol, vars = [22], markersize=1, label = "real current gfm+nn surrogate")
scatter!(p4_log, sol, vars = [23], markersize=1, label = "imag current gfm+nn surrogate")
p5 = plot(p1,p2,p1_log,p2_log,p3,p4,p3_log,p4_log, layout = (4,2), size = (1000,1000))
display_plots && display(p5)
##
################################# TRAINING #########################################
u₀ = x₀_nn  #stays same for a full training disturbance. 

function predict_gfm_nn(θ) 
    p = vcat(θ, p_inv, refs, p_fixed,  Vr0, Vi0) 
    _prob = remake(gfm_nn_prob, p=p, u0=u₀)    
    sol2 = solve(_prob, solver,  abstol=abstol, reltol=reltol, saveat=tsteps_train, 
                save_idxs=[i__ir_out, i__ii_out, i__ir_filter, i__ii_filter, i__ir_nn, i__ii_nn], #first two for loss function, rest for plotting
                 sensealg = ForwardDiffSensitivity())  
    return Array(sol2)
end
  
function loss_gfm_nn(θ)
    pred = predict_gfm_nn(θ)
    #To introduce batching, randomly select indices to calculate loss for
    loss = (mae(pred[1,:], ode_data_train[1,:]) / Ir_scale)/2 +
           (mae(pred[2,:], ode_data_train[2,:]) / Ii_scale)/2  
    return loss, pred, θ
end

function loss_fromdata(real_solution, predicted_solution)
    loss = (mae(predicted_solution[1,:], real_solution[1,:]) / Ir_scale)/2 + 
           (mae(predicted_solution[2,:], real_solution[2,:]) / Ii_scale)/2  
end 
    

list_plots = []
list_losses = Float64[]
list_θ = []
surr_data = []
list_gradnorm = Float64[]
iteration = 0

#Callback run extra time at end if you reach maxiters 
cb_gfm_nn = function(p, l, pred, θ) 
    #DISPLAY LOSS AND PLOT
    global iteration += 1 
    grad_norm = Statistics.norm(ForwardDiff.gradient(x -> first(loss_gfm_nn(x)),θ), 2) #Better to have a training infrastructure that saves and passes gradient instead of re-calculating 
    push!(list_gradnorm, grad_norm)
    push!(list_losses,l)
    push!(list_θ, θ[1:end]) #need[1:end], not sure why
    println("iteration:  ", iteration, "  loss: ", l)
    cb_gfm_nn_plot(pred)

    #UPDATE REFERENCES AND INITIAL CONDITIONS - if we don't do the NL solve, might not need any of this? 
#=     x₀, refs_int, Vr0_int, Vi0_int = initialize_sys!(sys_init, "gen1") #TODO don't need this each time
    global refs = refs_int
    global Vr0 = Vr0_int
    global Vi0 = Vi0_int
    p = vcat(θ, p_inv, refs, p_fixed, Vr0, Vi0)  
    x₀_nn = vcat(x₀, 0.0, 0.0,  x₀[5], x₀[19])
    f = get_init_gfm_nn(p, x₀[5], x₀[19])
    global res = nlsolve(f, x₀_nn)
    @assert converged(res)
    global u₀ = res.zero =#

    (l > lb_loss) && return false  
    return true 
end

#= if is_restart  #need to deal with re-starting in the correct range. 
    @load "checkpoint/mymodel.bson" p_nn opt list_loss list_grad iter
    iter += 1
end =#

ranges = extending_ranges(steps, group_size)

list_losses_best = []
list_θ_best = []
p_start = p_nn
for range in ranges
    global ode_data_train = ode_data[:,range]
    global tsteps_train = tsteps[range]
    res_gfm = DiffEqFlux.sciml_train(loss_gfm_nn, p_start, optimizer, GalacticOptim.AutoForwardDiff(), cb=cb_gfm_nn, maxiters=maxiters )
    global p_start = res_gfm.u  
    push!(list_θ_best, p_start) 
    push!(list_losses_best, res_gfm.minimum) 
    println("finished with range:", range)
    global iteration = 0 
end 

println("avg model achieves loss of ", loss_fromdata(ode_data, avgmodel_data))
println("surr model achieves loss of ", loss_gfm_nn(list_θ[end])[1]) #if maxiters is hit, last entry might not be lowest loss
surr_data = predict_gfm_nn(list_θ[end])
println("surr model achieves loss of ", loss_fromdata(ode_data, surr_data))
pcomp = plot_compare(ode_data, avgmodel_data, surr_data)
display(pcomp)

############################ SAVING PLOTS AND DATA ####################################
anim = Animation()
for plt in list_plots
    frame(anim, plt)
end

gif(anim, string("figs/", label ,"_train.gif"), fps = 100)
png(plot(list_losses), string("figs/",label,"_loss.png"))
writedlm( string("figs/", label ,"_loss.txt"), list_losses, ',')
png(plot(list_gradnorm), string("figs/", label,"_gradnorm.png")) 
png(list_plots[end], string("figs/", label ,"_final.png"))
png(pcomp, string("figs/", label ,"_comp.png"))
writedlm( string("figs/", label ,"_comp.txt"),  [tsteps'; ode_data; surr_data;avgmodel_data]', ',')
