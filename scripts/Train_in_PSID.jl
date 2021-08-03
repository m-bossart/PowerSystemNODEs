using Pkg
Pkg.activate(".")
using Revise
using Distributions
using OrdinaryDiffEq
using PowerSystems
using PowerSimulationsDynamics
const PSID = PowerSimulationsDynamics
const PSY = PowerSystems
using Plots
using FFTW
using Statistics
using NLsolve
include("../models/DynamicComponents.jl")
include("../models/InverterModels.jl")
include("../models/StaticComponents.jl")
include("../models/utils.jl")
include("../models/init_functions.jl")
########################PARAMETERS##############################################

#Make the non-reference bus a PV bus.
#set_bustype!(collect(get_components(Bus, sys_full, x->x.bustype==BusTypes.PQ ))[1],BusTypes.PV)

const train_split = 0.99    #proportion of faults for training (rest for test)
n_devices = 10             #number of devices in the surrogate
param_range = (0.9,1.1)
total_rating = 150.0  #MVA rating at the surrogate bus.

solver = Rodas5()
dtmax = 0.02
tspan = (0.0, 5.5)
step = 1e-2
tsteps = tspan[1]:step:tspan[2]
################################################################################
#NAMING FOR THE VARIOUS PSID SYSTEMS USED IN THIS SCRIPT
    #sys_faults: system with only Sources and PeriodicVariableSources for various faults
    #sys_surrogate: THe strucutre of the surrogate system... Do we need this?
    #sys_full: The structure of the full system
    #sys_train: Train system for the full system
    #sys_test: Test system for the full system.


###############BUILD THE SURROGATE SYSTEMS (FOR INITIALIZATION)#################
#PQbus = [b for b in get_components(Bus, sys_surrogate) if b.bustype == BusTypes.PQ][1]
#g = StaticGen("1", PQbus)
#set_rating!(g, total_rating/get_base_power(g))
#add_component!(sys_surrogate, g)
#inv_typ = inv_case78(get_name(g))
#add_component!(sys_surrogate, inv_typ, g)
#(sys_train_surrogate, sys_test_full) = build_train_test(sys_faults, sys_surrogate, train_split)

################BUILD THE FULL SYSTEMS (FOR GENERATING TRUTH DATA)##############
sys_faults = System("systems/fault_library.json")
sys_full = System("systems/base_system.json")
sys_train, sys_test = build_train_test(sys_faults, sys_full, 2, train_split, add_pvs = false) #Don't add PVS because can't build a sim with it yet.
@info "training set size:", length(collect(get_components(Source,sys_train)))
@info "test set size:", length(collect(get_components(Source,sys_test)))

to_json(sys_train,"systems/sys_train.json", force = true )
to_json(sys_test,"systems/sys_test.json", force = true)

##############################TRAINING##########################################
available_source = activate_next_source!(sys_train)

#Bus voltage is used in power flow, not source voltage. Need to set bus voltage from soure internal voltage
Vsource = get_internal_voltage(available_source)
set_magnitude!(get_bus(available_source),Vsource)
θsource = get_internal_angle(available_source)
set_angle!(get_bus(available_source),θsource)

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

#TODO The Source has the available field... get the available source and then get the dynamic injector from it.
active_source = collect(get_components(Source, sys_train,  x -> PSY.get_available(x)))[1]
V, θ = Source_to_function_of_time(active_source)
M = MassMatrix(19, 0)
gfm_func = ODEFunction(gfm, mass_matrix = M)


p_all = vcat(p_inv, refs, get_x(transormer), get_r(transformer))
##

dx = similar(x₀)
gfm(dx,x₀,p_all,0.0)
@assert all(isapprox.(dx, 0.0; atol=5e-6)) #TODO add the transformer impedance into gfm formulation so it can match!
f = get_init_gfm(p_all, x₀[5], x₀[19])
res = nlsolve(f, x₀)
gfm(dx,res.zero,p_all,0.0)
@assert all(isapprox.(dx, 0.0; atol=1e-8))


#TODO use x0 as a starting guess to initialize the various surrogates!
#TODO confirm the surrogate matches in a steady state sim (cannot use PVS yet!)




#refs = get_ext(collect(get_components(DynamicInverter, sys_surrogate))[1])["control_refs"]
#TODO : references should not be learned. Need to be included in the function throug a closrue?
#p_all = [p_start, p_refs ]

gfm_prob = ODEProblem(gfm_func,x₀,Float32.(tspan),p_inv)
#TODO Build an ODE problem with the various surrogate models
#TODO Solve the system using V,θ from above.



#TODO Will need to implement a method similar to get_activepower_series() for the PVS.
    #Will require a compute_output_current method for the PVS.
    #I can implement this once Jose has implemented the PVS dynamic model.
#TODO Make sure the surrogate matches the PSID system when we have identical models...


## PSEUDO CODE BELOW
function predict_solution(θ)        #Start with the non-UODE case, just optimize parameters.
    u₀ = initialize_surrogate()
    Array(solve(prob_surrogate))
end

function loss(θ)
    pred = predict_solution(θ)
    loss = real_solution - pred
end

list_losses = Float64[]

#TODO Need to figure out how to modify funcitons/data from within the callback.
#Look for examples that do this sort of thing.
cb_iters = function (θ,l)
    push!(list_losses,l)    #record loss

    if (loss < ϵ) # or if you've reach some number of epochs?
        activate_next_source!(sys)
        sim = Simulation!(MassMatrixModel, sys, pwd(),tspan)
        execture!(sim)
        #TODO: Modify the ground truth data here
        active_pvs = get_components(PeriodicVariableSource, sys,  x -> PSY.get_available(x))
        V, θ = PVS_to_function_of_time(active_pvs)
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
    return false
end

res = @time DiffEqFlux.sciml_train(loss,  p_nn,ADAM(0.1), cb = cb_iters, maxiters = 10) #TODO make maxiters = (#oftrainingfaults x iters_per_fault)
res = @time DiffEqFlux.sciml_train(loss,  p_nn,ADAM(0.1), cb = cb_threshold, maxiters = 10 )
