#Generates data for comparing various types of physics-based surrogates 
using Revise
using PowerSystems
using PowerFlows
using PowerSimulationNODE
using PowerSimulationsDynamicsSurrogates
const PSIDS = PowerSimulationsDynamicsSurrogates
using Logging
using Serialization
using Plots
include(joinpath(@__DIR__, "..", "build_datasets", "utils.jl"))
include(joinpath(@__DIR__, "..", "hpc_train", "utils.jl"))
system_name = "36bus_fix"
project_folder = "PowerSystemNODEs"
scratch_path = joinpath(pwd(), "..")

function show_surrogate_output_MW(p)
    @show "MW load", p[5] * p[1]
    @show "MW gfl", p[12] * p[2]    
    @show "MW gfm", p[33] * p[3]   
    @show "p_net", p[5] * p[1] + p[12] * p[2]+ p[33] * p[3]
end 

function halve_load(p)
    p_new = copy(p)
    MW_load = p[5] * p[1]
    MW_add = MW_load/4
    gfl_add = MW_add /p[12]
    gfm_add = MW_add /p[33]
    p_new[1] = p[1]/2    
    p_new[4] = p[4]/2    
    p_new[2] = p[2] + gfl_add 
    p_new[3] = p[3] + gfm_add 
    return p_new
end 

function remove_load(p)
    p_new = copy(p)
    MW_load = p[5] * p[1]
    MW_add = (MW_load * 0.99)/2     
    gfl_add = MW_add /p[12]
    gfm_add = MW_add /p[33]
    p_new[1] = p[1] * 0.01    #fails if zero (singular)
    p_new[4] = p[4] * 0.01    #fails if zero (singular)
    p_new[2] = p[2] + gfl_add 
    p_new[3] = p[3] + gfm_add 
    return p_new
end 

function remove_gfl(p)
    p_new = copy(p)
    MW_gfl = p[12] * p[2]
    MW_add = MW_gfl/2
    load_add = MW_add /p[5]
    gfm_add = MW_add /p[33]
    p_new[2] = 0.0 
    p_new[1] = p[1] + load_add 
    p_new[3] = p[3] + gfm_add 
    return p_new
end 

function remove_gfm(p)
    p_new = copy(p)
    MW_gfm = p[33] * p[3]
    MW_add = MW_gfm/2
    load_add = MW_add /p[5]
    gfl_add = MW_add /p[12]
    p_new[3] = 0.0 
    p_new[1] = p[1] + load_add 
    p_new[2] = p[2] + gfl_add 
    return p_new
end 

function load_only(p)
    p_new = copy(p)
    p_net_MW = p[5] * p[1] + p[12] * p[2]+ p[33] * p[3]
    p_new[1] = p_net_MW / p[5]
    p_new[2] = 0.0
    p_new[3] = 0.0
    return p_new
end 

function gfl_only(p)
    p_new = copy(p)
    p_net_MW = p[5] * p[1] + p[12] * p[2]+ p[33] * p[3]
    p_new[1] = 0.00001
    p_new[5] = 0.00001
    p_new[2] = p_net_MW / p[12]
    p_new[3] = 0.0
    return p_new
end 
function gfm_only(p)
    p_new = copy(p)
    p_net_MW = p[5] * p[1] + p[12] * p[2]+ p[33] * p[3]
    p_new[1] = 0.00001
    p_new[5] = 0.00001
    p_new[2] = 0.0
    p_new[3] = p_net_MW / p[33]
    return p_new
end 

function modify_trainable(p, percent)
    p_new = copy(p)
    p_new[1:4] .= p[1:4] .* (1+percent)     
    p_new[15:16] .= p[15:16]  .* (1+percent)
    p_new[18:19] .= p[18:19]  .* (1+percent)
    p_new[36] = p[36] * (1+percent)
    p_new[38] = p[38] * (1+percent)
end 

if Sys.iswindows() || Sys.isapple()
    include(joinpath(pwd(), "scripts", "hpc_train", "utils.jl"))
    sys = System(joinpath("systems", string(system_name, ".json")))
    surrogate_buses = vcat(21:29, 21:39)
    #θ = determine_p_start(sys, surrogate_buses)
    #Serialization.serialize(joinpath("starting_parameters", "p_start_physics"), θ)
end

starting_parameter_options = []
experiment_names = [] 
θ = determine_p_start(sys, surrogate_buses)
push!(experiment_names, "exp_09_08_24_Load+GFL+GFM")
push!(starting_parameter_options, θ)
push!(experiment_names, "exp_09_08_24_Load+GFM")
push!(starting_parameter_options, remove_gfl(θ))
push!(experiment_names, "exp_09_08_24_Load+GFL")
push!(starting_parameter_options, remove_gfm(θ))
push!(experiment_names, "exp_09_08_24_GFL")
push!(starting_parameter_options, gfl_only(θ))
push!(experiment_names, "exp_09_08_24_GFM")
push!(starting_parameter_options, gfm_only(θ))

for (exp_name, starting_parameters) in zip(experiment_names, starting_parameter_options)
    train_folder = joinpath("data_from_hpc", exp_name)
    system_name = "36bus_fix"
    project_folder = "PowerSystemNODEs"
    scratch_path = joinpath(pwd(), "..")
    _copy_full_system_to_train_directory(
        scratch_path,
        project_folder,
        train_folder,
        system_name,
    )
    p = TrainParams(
        train_id = "BASE",
        surrogate_buses = vcat(21:29, 31:39),
        train_data = (
            id = "1",
            operating_points = repeat(
                [
                    RandomOperatingPointXiao(
                        generator_voltage_range = (0.94, 1.06),
                        generator_power_range = (0.0, 1.0),
                        load_multiplier_range = (0.5, 1.5),
                    ),
                ],
                1,
            ),
            perturbations = repeat(
                [[PSIDS.RandomLoadChange(time = 1.0, load_multiplier_range = (0.0, 2.0))]],
                1,
            ),
            params = PSIDS.GenerateDataParams(
                solver = "Rodas5",
                solver_tols = (reltol = 1e-3, abstol = 1e-6),
                tspan = (0.0, 10.0),
                tstops = 0.0:0.1:10.0,
                tsave = 0.0:0.1:10.0,
                formulation = "MassMatrix",
                all_branches_dynamic = false,
                all_lines_dynamic = false,
                seed = 11,
            ),
            system = "full",
        ),
        validation_data = (
            id = "1",
            operating_points = repeat(
                [
                    RandomOperatingPointXiao(
                        generator_voltage_range = (0.94, 1.06),
                        generator_power_range = (0.0, 1.0),
                        load_multiplier_range = (0.5, 1.5),
                    ),
                ],
                1,
            ),
            perturbations = repeat(
                [[PSIDS.RandomLoadChange(time = 1.0, load_multiplier_range = (0.0, 2.0))]],
                1,
            ),
            params = PSIDS.GenerateDataParams(
                solver = "Rodas5",
                solver_tols = (reltol = 1e-3, abstol = 1e-6),
                tspan = (0.0, 10.0),
                tstops = 0.0:0.1:10.0,
                tsave = 0.0:0.1:10.0,
                formulation = "MassMatrix",
                all_branches_dynamic = false,
                all_lines_dynamic = false,
                seed = 22,
            ),
        ),
        test_data = (
            id = "1",
            operating_points = repeat(
                [
                    RandomOperatingPointXiao(
                        generator_voltage_range = (0.94, 1.06),
                        generator_power_range = (0.0, 1.0),
                        load_multiplier_range = (0.5, 1.5),
                    ),
                ],
                100,
            ),
            perturbations = repeat(
                [[PSIDS.RandomLoadChange(time = 1.0, load_multiplier_range = (0.0, 2.0))]],
                1,
            ),
            params = PSIDS.GenerateDataParams(
                solver = "Rodas5",
                solver_tols = (reltol = 1e-3, abstol = 1e-6),
                tspan = (0.0, 10.0),
                tstops = 0.0:0.1:10.0,
                tsave = 0.0:0.1:10.0,
                formulation = "MassMatrix",
                all_branches_dynamic = false,
                all_lines_dynamic = false,
                seed = 33,
            ),
        ),
        model_params = MultiDeviceParams(name = "source_1"),
        optimizer = [
            (
                auto_sensealg = "ForwardDiff",
                algorithm = "Adam",
                log_η = -7.0,
                initial_stepnorm = 0.0,
                maxiters = 1,
                steadystate_solver = (
                    solver = "NLSolveJL",
                    reltol = 1e-4,
                    abstol = 1e-4,
                    termination = "RelSafeBest",
                ),
                dynamic_solver = (
                    solver = "Rodas5",
                    sensealg = "QuadratureAdjoint",
                    reltol = 1e-3,
                    abstol = 1e-6,
                    maxiters = 1e5,
                    force_tstops = true,
                ),
                lb_loss = 0.0,
                curriculum = "individual faults",
                curriculum_timespans = [(tspan = (0.0, 10.0), batching_sample_factor = 1.0)],
                fix_params = [   #fix zip to constant impedance 
                    #:P_fraction_1,
                    #:P_fraction_2,
                    #:P_fraction_3,
                    #:Q_fraction_1,
                    :base_power_zip,
                    :max_active_power_Z,
                    :max_active_power_I,
                    :max_active_power_P,
                    :max_reactive_power_Z,
                    :max_reactive_power_I,
                    :max_reactive_power_P,
                    :base_power_gfl,
                    :rated_voltage_gfl,
                    :rated_current_gfl,
                    #:Kp_p,
                    #:Ki_p,
                    :ωz_gfl,
                    #:Kp_q,
                    #:Ki_q,
                    :ωf_gfl,
                    :kpc_gfl,
                    :kic_gfl,
                    :kffv_gfl,
                    :voltage_gfl,
                    :ω_lp,
                    :kp_pll,
                    :ki_pll,
                    :lf_gfl,
                    :rf_gfl,
                    :cf_gfl,
                    :lg_gfl,
                    :rg_gfl,
                    :base_power_gfm,
                    :rated_voltage_gfm,
                    :rated_current_gfm,
                    #:Rp,    
                    :ωz_gfm,
                    #:kq,
                    :ωf_gfm,
                    :kpv,
                    :kiv,
                    :kffv_gfm,
                    :rv,
                    :lv,
                    :kpc_gfm,
                    :kic_gfm,
                    :kffi,
                    :ωad,
                    :kad,
                    :voltage_gfm,
                    :lf_gfm,
                    :rf_gfm,
                    :cf_gfm,
                    :lg_gfm,
                    :rg_gfm,
                ],
                loss_function = (α = 0.5, β = 1.0, residual_penalty = 1.0e2),
            ),
        ],
        check_validation_loss_iterations = [10],
        p_start = starting_parameters,
        final_validation_loss = true,
        time_limit_buffer_seconds = 7200,
        rng_seed = 11,
        output_mode_skip = 1,
        train_time_limit_seconds = 1e9,
        base_path = joinpath(scratch_path, project_folder, train_folder),
        system_path = joinpath(
            scratch_path,
            project_folder,
            train_folder,
            PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
            string(system_name, ".json"),
        ),
    )
    build_subsystems(p)
    mkpath(joinpath(p.base_path, PowerSimulationNODE.INPUT_FOLDER_NAME))

    if exp_name == "exp_09_08_24_Load+GFL+GFM"
        generate_train_data(p)
        generate_validation_data(p)
        generate_test_data(p)       #Can just copy over from a prior experiment; this takes a while to run...
    else 
        cp(
            joinpath(pwd(), "data_from_hpc", "exp_09_08_24_Load+GFL+GFM", "input_data", "train_data_1"),
            p.train_data_path
        )
        cp(
            joinpath(pwd(), "data_from_hpc", "exp_09_08_24_Load+GFL+GFM", "input_data", "validation_data_1"),
            p.validation_data_path
        )
        cp(
            joinpath(pwd(), "data_from_hpc", "exp_09_08_24_Load+GFL+GFM", "input_data", "test_data_1"),
            p.test_data_path
        )
    end 
    train(p)
    input_param_file = joinpath(p.base_path, "input_data", "train_001.json")
    PowerSimulationNODE.serialize(p, input_param_file)
end  

##
##########################
# Visualize datasets (use before attempting to train)
#train_dataset = Serialization.deserialize(p.train_data_path)
#display(visualize_dataset(train_dataset))
#validation_dataset = Serialization.deserialize(p.validation_data_path)
#display(visualize_dataset(validation_dataset))
#test_dataset = Serialization.deserialize(p.test_data_path)
#display(visualize_dataset(test_dataset))
#
#sys_validation = PowerSimulationNODE.node_load_system(p.surrogate_system_path)
#sys_validation_aux = deepcopy(sys_validation)
#PowerSimulationNODE.add_surrogate_psid!(sys_validation, p.model_params, train_dataset)
#to_json(sys_validation, p.modified_surrogate_system_path, force = true)
#
#
#
#sys = System(p.system_path)
#p_start =  determine_p_start(sys, p.surrogate_buses)
#show_surrogate_output_MW(p_start)
#p_start = gfm_only(p_start)
##p_start = remove_load(p_start)
#show_surrogate_output_MW(p_start)
#
##show_surrogate_output_MW(p_without_load)
##p_start = remove_gfl(p_start)
##show_surrogate_output_MW(p_without_gfl)
##p_start = remove_gfm(p_start)
##show_surrogate_output_MW(p_without_gfm)
#
#p_start[5] * p_start[1] + p_start[12] * p_start[2]+ p_start[33] * p_start[3]
#
#@assert isapprox(p_start[5] * p_start[1] + p_start[12] * p_start[2]+ p_start[33] * p_start[3], -92.73999999999987)
#validation_sys = PowerSimulationNODE.node_load_system(p.modified_surrogate_system_path)
#validation_sys_aux = PowerSimulationNODE.node_load_system(p.surrogate_system_path)
#data_collection_location = Serialization.deserialize(p.data_collection_location_path)[2]
#
#surrogate_dataset = generate_surrogate_dataset(
#    validation_sys,
#    validation_sys_aux,
#    p_start,
#    validation_dataset,
#    p.validation_data,
#    data_collection_location,
#    p.model_params,
#)
#
#@show validation_dataset[1].device_terminal_data["Bus 6 -> Bus 26"][:ir][1]
#@show surrogate_dataset[1].device_terminal_data["Bus 6 -> Bus 26"][:ir][1]
#@show validation_dataset[1].device_terminal_data["Bus 6 -> Bus 26"][:ii][1]
#@show surrogate_dataset[1].device_terminal_data["Bus 6 -> Bus 26"][:ii][1]
##
######################################################################################
####################################### TRAIN ########################################
######################################################################################

######################################################################################
############################### ANALYZE AND VISUALIZE ################################
######################################################################################
##

#visualize_training(input_param_file, collect(1:20))
##
#animate_training(input_param_file, skip = 1)
#a = generate_summary(joinpath(p.base_path, "output_data"))
#pp = visualize_summary(a)
#print_high_level_output_overview(a, p.base_path)
