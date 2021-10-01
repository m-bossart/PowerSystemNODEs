include("../models/constants.jl")
include("../models/DynamicComponents.jl")
include("../models/InverterModels.jl")
include("../models/utils.jl")
include("../models/parameter_utils.jl")
include("../models/init_functions.jl")
configure_logging(console_level = Logging.Error)

################BUILD THE TRAINING SYSTEMS FOR GENERATING TRUTH DATA#############
sys_faults = System("systems/fault_library_3invs_vsms_20%lossP.json")
sys_full = System("systems/base_system_3invs_vsms_20%lossP.json")
sys_train, sys_test = build_train_test(sys_faults, sys_full, 2, train_split, add_pvs = true)
@info "training set size:", length(collect(get_components(Source,sys_train)))
@info "test set size:", length(collect(get_components(Source,sys_test)))
to_json(sys_train,"systems/sys_train.json", force = true)
to_json(sys_test,"systems/sys_test.json", force = true)


############################# GENERATE TRUE SOLUTION ###########################
available_source = activate_next_source!(sys_train)
set_bus_from_source(available_source) #Bus voltage is used in power flow, not source voltage. Need to set bus voltage from soure internal voltage

sim_full = Simulation!(
    MassMatrixModel,
    sys_train,
    pwd(),
    tspan,
)

@info "train system power flow", solve_powerflow(sys_train)["flow_results"]
@info "train system power flow", solve_powerflow(sys_train)["bus_results"]
show_states_initial_value(sim_full)

println("time to solve train system for generating truth data:")
@time execute!(sim_full,
        solver,
        abstol = abstol,
        reltol = reltol,
        reset_simulation=false, saveat=tsteps );
#tsteps = sim.solution.t  #to use the solver time steps, uncomment this line and get rid of saveat above 

sol_full = read_results(sim_full)

active_source = collect(get_components(Source, sys_train,  x -> PSY.get_available(x)))[1]
Vmag_bus = get_voltage_magnitude_series(sol_full, 2)
θ_bus = get_voltage_angle_series(sol_full, 2)
Vmag_internal = get_state_series(sol_full, ("source1",:Vt))
θ_internal = get_state_series(sol_full, ("source1",:θt))

Vr_scale = 1/(maximum(Vmag_bus[2].*cos.(θ_bus[2])) - minimum(Vmag_bus[2].*cos.(θ_bus[2])))
Vi_scale = 1/(maximum(Vmag_bus[2].*sin.(θ_bus[2])) - minimum(Vmag_bus[2].*sin.(θ_bus[2])))

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
plot(p1,p2,p1_log,p2_log, layout=(2,2), size = (1000,1000))

#plot the other bus voltage too? to see if it is significantly different?? 
ode_data = get_total_current_series(sim_full) #TODO Better to measure current at the PVS (implement method after PVS is complete)
p3 = plot(tsteps, ode_data[1,:], label = "Ir full model - PSID")
p4 = plot(tsteps, ode_data[2,:], label = "Ii full model - PSID")
p3_log = plot(tsteps,ode_data[1,:], label =  "Ir full model - PSID", xaxis=:log)
p4_log = plot(tsteps,ode_data[2,:], label =  "Ii full model - PSID", xaxis=:log)

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
sim_simp = Simulation!(
    MassMatrixModel,
    sys_init,
    pwd(),
    tspan,
)
@info "initialize system power flow", solve_powerflow(sys_init)["flow_results"]
@info "initialize system power flow", solve_powerflow(sys_init)["bus_results"]
show_states_initial_value(sim_simp)
@time execute!(sim_simp,
        solver,
        abstol = abstol,
        reltol = reltol,
        reset_simulation=false, saveat=tsteps);


avgmodel_data_p = get_real_current_series(read_results(sim_simp),"gen1")
avgmodel_data = get_total_current_series(sim_simp)

plot!(p3, tsteps, avgmodel_data[1,:], markersize=1, label = "Ir simple model - PSID")
plot!(p4, tsteps, avgmodel_data[2,:], markersize=1, label = "Ii simple model - PSID")
plot!(p3_log, tsteps, avgmodel_data[1,:], markersize=1, label = "Ir simple model - PSID")
plot!(p4_log, tsteps, avgmodel_data[2,:], markersize=1, label = "Ii simple model - PSID")


##### INITIALIZE THE GFM+NN SURROGATE AND BUILD THE TRAINING PROBLEM ###########
nn_org = nn_scale 
nn_scale = 0

nn = build_nn(4, 2, nn_width, nn_hidden, nn_activation) #include currents 
p_nn = initial_params(nn)
n_weights_nn = length(p_nn)
p_all = vcat(p_nn, p_inv, refs, p_fixed, Vr0, Vi0)
x₀_nn = vcat(x₀, 0.0, 0.0, x₀[5], x₀[19])

h = get_init_gfm_nn(p_all, x₀[5], x₀[19])
res_nn= nlsolve(h, x₀_nn)
@assert converged(res_nn)
dx = similar(x₀_nn)
gfm_nn(dx,res_nn.zero,p_all,0.0)
@assert all(isapprox.(dx, 0.0; atol=1e-8))

M = MassMatrix(21, 2)
gfm_nn_func = ODEFunction(gfm_nn, mass_matrix = M)
gfm_nn_prob = ODEProblem(gfm_nn_func, x₀_nn, tspan, p_all)
sol = solve(gfm_nn_prob, solver, abstol=abstol, reltol=reltol,  saveat=tsteps )
plot!(p3, sol.t,sol[22,:], label = "Ir simple model - hand")
plot!(p4, sol.t, sol[23,:], label = "Ii simple model - hand")

nn_scale = nn_org

#UNCOMMENT TO INCLUDE THE UNTRAINED SURROGATE IN THE COMPARISON
sol = solve(gfm_nn_prob, solver, abstol=abstol, reltol=reltol,  saveat=tsteps )
#plot!(p3, sol.t,sol[22,:], label = "Ir surrogate (untrained)")
#plot!(p4, sol.t, sol[23,:], label = "Ii surrogate (untrained)")

if (plot_log)
    p5 = plot(p1,p2,p1_log,p2_log,p3,p4,p3_log,p4_log, layout = (4,2))
else 
    p5 = plot(p1,p2,p3,p4, layout = (2,2))
end
display_plots && display(p5)

##
################################# TRAINING #########################################
u₀ = x₀_nn  #stays same for a full training disturbance. 
batch = range(1,length = group_size)

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
    @show batch
    if(loss_function == "mae")
        loss = (mae(pred[1,batch], ode_data_train[1,batch]) / Ir_scale)/2 +
            (mae(pred[2,batch], ode_data_train[2,batch]) / Ii_scale)/2  

    elseif(loss_function == "mse")
        loss = (mse(pred[1,batch], ode_data_train[1,batch]) / Ir_scale)/2 +
        (mse(pred[2,batch], ode_data_train[2,batch]) / Ii_scale)/2  
    end
    return loss, pred, θ
end

function loss_fromdata(real_solution, predicted_solution)
    if(loss_function == "mae")
        loss = (mae(predicted_solution[1,:], real_solution[1,:]) / Ir_scale)/2 + 
            (mae(predicted_solution[2,:], real_solution[2,:]) / Ii_scale)/2  
    elseif(loss_function == "mse")
        loss = (mse(predicted_solution[1,:], real_solution[1,:]) / Ir_scale)/2 + 
            (mse(predicted_solution[2,:], real_solution[2,:]) / Ii_scale)/2 
    end 
end 
    

list_plots = []
list_losses = Float64[]
list_θ = []
surr_data = []
#list_gradnorm = Float64[]
iteration = 0

#Callback run extra time at end if you reach maxiters 
cb_gfm_nn = function(p, l, pred, θ) 
    #DISPLAY LOSS AND PLOT
    global iteration += 1 
    #grad_norm = Statistics.norm(ForwardDiff.gradient(x -> first(loss_gfm_nn(x)),θ), 2) #Better to have a training infrastructure that saves and passes gradient instead of re-calculating 
    #push!(list_gradnorm, grad_norm)
    push!(list_losses,l)
    push!(list_θ, θ[1:end]) #need[1:end], not sure why
    println("iteration:  ", iteration, "  loss: ", l)
    cb_gfm_nn_plot(pred, plot_log)

    #UPDATE REFERENCES AND INITIAL CONDITIONS - 
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
    length_curr = length(ode_data_train[1,:])
    if (batching_factor == 1)
        global batch = range(1,length=length_curr)
    else
        global batch = rand(range(1,length=length_curr),Int(length_curr/batching_factor))
    end 
    (l > lb_loss) && return false  
    return true 
end


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

#Use the full range for reporting results 
ode_data_train = ode_data
tsteps_train = tsteps
batch = range(1,length=length(ode_data_train[1,:]))

@show batch

println("avg model achieves loss of ", loss_fromdata(ode_data, avgmodel_data))
#if maxiters is hit, last entry might not be lowest loss
#surr_data = predict_gfm_nn(list_θ[argmin(list_losses)])
sublist_θ = list_θ[end-maxiters+1:end-1]
sublist_losses = list_losses[end-maxiters+1:end-1] 
best_θ = sublist_θ[argmin(sublist_losses)]
surr_data = predict_gfm_nn(best_θ)
println("surr model achieves loss of ", loss_gfm_nn(best_θ)[1])
pcomp = plot_compare(ode_data, avgmodel_data, surr_data)
display(pcomp)

############################ SAVING PLOTS AND DATA ####################################

png(plot(list_losses), string("figs/",label,"_loss.png"))
writedlm( string("figs/", label ,"_loss.txt"), list_losses, ',')
#png(plot(list_gradnorm), string("figs/", label,"_gradnorm.png")) 
png(list_plots[end], string("figs/", label ,"_final.png"))
png(pcomp, string("figs/", label ,"_comp.png"))
writedlm( string("figs/", label ,"_comp.txt"),  [tsteps'; ode_data; surr_data;avgmodel_data]', ',')
#SAVE BEST θ!
@save string("figs/",label) best_θ

anim = Animation()
for plt in list_plots
    frame(anim, plt)
end

gif(anim, string("figs/", label ,"_train.gif"), fps = 100)
