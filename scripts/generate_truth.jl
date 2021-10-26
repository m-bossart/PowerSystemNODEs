#optimizer = ADAM(0.01)# BFGS(initial_stepnorm=0.01)# ADAM(0.01)  #BFGS(), Optim.KrylovTrustRegion(),
train_split = 0.9999
solver = Rodas4()       #KenCarp4(), QBDF(), TRBDF2() 
abstol = 1e-6
reltol = 1e-3
tfault = 0.1
tspan = (0.0, 2.0)
steps = 150
group_size = 150
tsteps = 10 .^ (range(log10(tfault), log10(tspan[2]), length = steps))
tsteps = tspan[1]:((tspan[2] - tspan[1]) / steps):tspan[2]
lb_loss = 0
nn_width = 5
nn_hidden = 5
maxiters = 5#250  #Change 
nn_activation = gelu
nn_scale = 1.0
plot_log = false
display_plots = false
loss_function = "mae"

################BUILD THE TRAINING SYSTEMS FOR GENERATING TRUTH DATA#############
sys_faults = System("systems/fault_library_3invs_vsms_20%lossP.json")
sys_full = System("systems/base_system_3invs_vsms_20%lossP.json")
sys_train, sys_test = build_train_test(sys_faults, sys_full, 2, train_split, add_pvs = true)
@info "training set size:", length(collect(get_components(Source, sys_train)))
@info "test set size:", length(collect(get_components(Source, sys_test)))
to_json(sys_train, "systems/sys_train.json", force = true)
to_json(sys_train, "input_data/system.json", force = true)
to_json(sys_test, "systems/sys_test.json", force = true)

############################# GENERATE TRUE SOLUTION ###########################
available_source = activate_next_source!(sys_train)
set_bus_from_source(available_source) #Bus voltage is used in power flow, not source voltage. Need to set bus voltage from soure internal voltage

sim_full = Simulation!(MassMatrixModel, sys_train, pwd(), tspan)
res = small_signal_analysis(sim_full)
Stiffness =
    maximum(abs.(real(res.eigenvalues[1:(end - 2)]))) /
    minimum(abs.(real(res.eigenvalues[1:(end - 2)])))
@info "Stiffness", Stiffness

@info "train system power flow", solve_powerflow(sys_train)["flow_results"]
@info "train system power flow", solve_powerflow(sys_train)["bus_results"]
show_states_initial_value(sim_full)

println("time to solve train system for generating truth data:")
@time execute!(
    sim_full,
    solver,
    abstol = abstol,
    reltol = reltol,
    reset_simulation = false,
    saveat = tsteps,
);
#tsteps = sim.solution.t  #to use the solver time steps, uncomment this line and get rid of saveat above 

sol_full = read_results(sim_full)

active_source = collect(get_components(Source, sys_train, x -> PSY.get_available(x)))[1]
Vmag_bus = get_voltage_magnitude_series(sol_full, 2)
θ_bus = get_voltage_angle_series(sol_full, 2)
Vmag_internal = get_state_series(sol_full, ("source1", :Vt))
θ_internal = get_state_series(sol_full, ("source1", :θt))

Vr_scale =
    1 / (maximum(Vmag_bus[2] .* cos.(θ_bus[2])) - minimum(Vmag_bus[2] .* cos.(θ_bus[2])))
Vi_scale =
    1 / (maximum(Vmag_bus[2] .* sin.(θ_bus[2])) - minimum(Vmag_bus[2] .* sin.(θ_bus[2])))

ode_data = get_total_current_series(sim_full) #TODO Better to measure current at the PVS (implement method after PVS is complete)

#################### BUILD INITIALIZATION SYSTEM ###############################
sys_init, p_inv = build_sys_init(sys_train) #returns p_inv, the set of average parameters 
transformer = collect(get_components(Transformer2W, sys_init))[1]
pvs = collect(get_components(PeriodicVariableSource, sys_init))[1]
p_fixed = [get_x(transformer) + get_X_th(pvs), get_r(transformer) + get_R_th(pvs)]
x₀, refs, Vr0, Vi0 = initialize_sys!(sys_init, "gen1")
Vm, Vθ = Source_to_function_of_time(get_dynamic_injector(active_source))
p_ode = vcat(p_inv, refs, p_fixed)
sim_simp = Simulation!(MassMatrixModel, sys_init, pwd(), tspan)
@info "initialize system power flow", solve_powerflow(sys_init)["flow_results"]
@info "initialize system power flow", solve_powerflow(sys_init)["bus_results"]
show_states_initial_value(sim_simp)
@time execute!(
    sim_simp,
    solver,
    abstol = abstol,
    reltol = reltol,
    reset_simulation = false,
    saveat = tsteps,
);

avgmodel_data_p = get_real_current_series(read_results(sim_simp), "gen1")
avgmodel_data = get_total_current_series(sim_simp)

#build dataframe: #Change to dictionary, store as json instead of arrow. 
#= df = DataFrame(
    :id => [get_name(pvs)],
    :tsteps => [tsteps],
    :ir_true => [ode_data[1, :]],
    :ii_true => [ode_data[2, :]],
    :ir_ver => [avgmodel_data[1, :]],      
    :ii_ver => [avgmodel_data[2, :]],
    :p_ode => [p_ode],
    :x₀ => [x₀],
    :V₀ => [[Vr0, Vi0]],
) =#

d = Dict{String, Dict{Symbol, Any}}()
d[get_name(pvs)] = Dict(
    :tsteps => tsteps,
    :ir_ground_truth => ode_data[1, :],
    :ii_ground_truth => ode_data[2, :],
    :ir_node_off => avgmodel_data[1, :],
    :ii_node_off => avgmodel_data[2, :],
    :p_ode => p_ode,
    :x₀ => x₀,
    :V₀ => [Vr0, Vi0],
)

open("input_data/data.json", "w") do io
    JSON3.write(io, d)
end

default_params = NODETrainParams()
open("train_parameters/default_NODE_params.json", "w") do io
    JSON3.write(io, default_params)
end
