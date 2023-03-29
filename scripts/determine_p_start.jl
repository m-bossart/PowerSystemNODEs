using PowerSystems
using PowerSimulationsDynamics
using PowerSimulationNODE
using PowerSimulationsDynamicsSurrogates
const PSIDS = PowerSimulationsDynamicsSurrogates
const PSN = PowerSimulationNODE
using Serialization
using Plots

include("build_datasets/utils.jl")
include("hpc_train/utils.jl")
train_folder = "exp_starting_parameters"
system_name = "36Bus_CR"
project_folder = "PowerSystemNODEs"
scratch_path = joinpath(pwd(), "..")

_copy_full_system_to_train_directory(
    scratch_path,
    project_folder,
    train_folder,
    system_name,
)

sys = System(joinpath("systems", string(system_name, ".json")))
surrogate_buses = [21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37, 38, 39]

p_powers = zeros(6)
p_gfl = zeros(20)
p_gfm = zeros(22)
p_load = zeros(6)

#CALCULATE AVERAGE GFL PARAMETERS 
n_gfl = 0
for i in get_components(
    DynamicInverter{
        AverageConverter,
        OuterControl{ActivePowerPI, ReactivePowerPI},
        CurrentModeControl,
        FixedDCSource,
        KauraPLL,
        LCLFilter,
    },
    sys,
)
    static_injector = get_component(StaticInjection, sys, get_name(i))
    if get_number(get_bus(static_injector)) in surrogate_buses
        println(get_name(static_injector))
        n_gfl += 1
        converter = get_converter(i)
        p_gfl[PSN.gfl_indices[:params][:rated_voltage_gfl]] += get_rated_voltage(converter)
        p_gfl[PSN.gfl_indices[:params][:rated_current_gfl]] += get_rated_current(converter)

        active_power = PSY.get_active_power_control(get_outer_control(i))
        p_gfl[PSN.gfl_indices[:params][:Kp_p]] += get_Kp_p(active_power)
        p_gfl[PSN.gfl_indices[:params][:Ki_p]] += get_Ki_p(active_power)
        p_gfl[PSN.gfl_indices[:params][:ωz_gfl]] += get_ωz(active_power)

        reactive_power = PSY.get_reactive_power_control(get_outer_control(i))
        p_gfl[PSN.gfl_indices[:params][:Kp_q]] += get_Kp_q(reactive_power)
        p_gfl[PSN.gfl_indices[:params][:Ki_q]] += get_Ki_q(reactive_power)
        p_gfl[PSN.gfl_indices[:params][:ωf_gfl]] += get_ωf(reactive_power)

        inner_control = get_inner_control(i)
        p_gfl[PSN.gfl_indices[:params][:kpc_gfl]] += get_kpc(inner_control)
        p_gfl[PSN.gfl_indices[:params][:kic_gfl]] += get_kic(inner_control)
        p_gfl[PSN.gfl_indices[:params][:kffv_gfl]] += get_kffv(inner_control)

        dc_source = get_dc_source(i)
        p_gfl[PSN.gfl_indices[:params][:voltage_gfl]] += get_voltage(dc_source)

        freq_estimator = get_freq_estimator(i)
        p_gfl[PSN.gfl_indices[:params][:ω_lp]] += get_ω_lp(freq_estimator)
        p_gfl[PSN.gfl_indices[:params][:kp_pll]] += get_kp_pll(freq_estimator)
        p_gfl[PSN.gfl_indices[:params][:ki_pll]] += get_ki_pll(freq_estimator)

        filter = get_filter(i)
        p_gfl[PSN.gfl_indices[:params][:lf_gfl]] += get_lf(filter)
        p_gfl[PSN.gfl_indices[:params][:rf_gfl]] += get_rf(filter)
        p_gfl[PSN.gfl_indices[:params][:cf_gfl]] += get_cf(filter)
        p_gfl[PSN.gfl_indices[:params][:lg_gfl]] += get_lg(filter)
        p_gfl[PSN.gfl_indices[:params][:rg_gfl]] += get_rg(filter)
    end
end
p_gfl = p_gfl ./ n_gfl

n_gfm = 0
for i in get_components(
    DynamicInverter{
        AverageConverter,
        OuterControl{ActivePowerDroop, ReactivePowerDroop},
        VoltageModeControl,
        FixedDCSource,
        FixedFrequency,
        LCLFilter,
    },
    sys,
)
    static_injector = get_component(StaticInjection, sys, get_name(i))
    if get_number(get_bus(static_injector)) in surrogate_buses
        println(get_name(static_injector))
        n_gfm += 1
        converter = get_converter(i)
        p_gfm[PSN.gfm_indices[:params][:rated_voltage_gfm]] += get_rated_voltage(converter)
        p_gfm[PSN.gfm_indices[:params][:rated_current_gfm]] += get_rated_current(converter)

        active_power = PSY.get_active_power_control(get_outer_control(i))
        p_gfm[PSN.gfm_indices[:params][:Rp]] += get_Rp(active_power)
        p_gfm[PSN.gfm_indices[:params][:ωz_gfm]] += get_ωz(active_power)

        reactive_power = PSY.get_reactive_power_control(get_outer_control(i))
        p_gfm[PSN.gfm_indices[:params][:kq]] += get_kq(reactive_power)
        p_gfm[PSN.gfm_indices[:params][:ωf_gfm]] += get_ωf(reactive_power)

        inner_control = get_inner_control(i)
        p_gfm[PSN.gfm_indices[:params][:kpv]] += get_kpv(inner_control)
        p_gfm[PSN.gfm_indices[:params][:kiv]] += get_kiv(inner_control)
        p_gfm[PSN.gfm_indices[:params][:kffv_gfm]] += get_kffv(inner_control)
        p_gfm[PSN.gfm_indices[:params][:rv]] += get_rv(inner_control)
        p_gfm[PSN.gfm_indices[:params][:lv]] += get_lv(inner_control)
        p_gfm[PSN.gfm_indices[:params][:kpc_gfm]] += get_kpc(inner_control)
        p_gfm[PSN.gfm_indices[:params][:kic_gfm]] += get_kpc(inner_control)
        p_gfm[PSN.gfm_indices[:params][:kffi]] += get_kffi(inner_control)
        p_gfm[PSN.gfm_indices[:params][:ωad]] += get_ωad(inner_control)
        p_gfm[PSN.gfm_indices[:params][:kad]] += get_kad(inner_control)

        dc_source = get_dc_source(i)
        p_gfm[PSN.gfm_indices[:params][:voltage_gfm]] += get_voltage(dc_source)

        freq_estimator = get_freq_estimator(i)

        filter = get_filter(i)
        p_gfm[PSN.gfm_indices[:params][:lf_gfm]] += get_lf(filter)
        p_gfm[PSN.gfm_indices[:params][:rf_gfm]] += get_rf(filter)
        p_gfm[PSN.gfm_indices[:params][:cf_gfm]] += get_cf(filter)
        p_gfm[PSN.gfm_indices[:params][:lg_gfm]] += get_lg(filter)
        p_gfm[PSN.gfm_indices[:params][:rg_gfm]] += get_rg(filter)
    end
end
p_gfm = p_gfm ./ n_gfm

n_load = 0
for i in get_components(StandardLoad, sys)
    if get_number(get_bus(i)) in surrogate_buses
        n_load += 1
        p_load[PSN.zip_indices[:params][:max_active_power_Z]] +=
            get_max_impedance_active_power(i)
        p_load[PSN.zip_indices[:params][:max_reactive_power_Z]] +=
            get_max_impedance_reactive_power(i)
        p_load[PSN.zip_indices[:params][:max_active_power_I]] +=
            get_max_current_active_power(i)
        p_load[PSN.zip_indices[:params][:max_reactive_power_I]] +=
            get_max_current_reactive_power(i)
        p_load[PSN.zip_indices[:params][:max_active_power_P]] +=
            get_max_constant_active_power(i)
        p_load[PSN.zip_indices[:params][:max_reactive_power_P]] +=
            get_max_constant_reactive_power(i)
    end
end
p_load = p_load ./ n_load
p_powers = [n_load, n_load, n_gfl, n_gfl, n_gfm, n_gfm] ./ (n_load + n_gfl + n_gfm) #hardcoded order...

p_start = vcat(p_powers, p_load, p_gfl, p_gfm)
##
p = TrainParams(
    train_id = "BASE",
    surrogate_buses = surrogate_buses,
    train_data = (
        id = "1",
        operating_points = PSIDS.SurrogateOperatingPoint[PSIDS.GenerationLoadScale(
            generation_scale = 1.0,
            load_scale = 1.0,
        ),],
        perturbations = repeat(
            [[PSIDS.RandomLoadChange(time = 1.0, load_multiplier_range = (0.0, 2.0))]],
            3,
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
            seed = 1,
        ),
        system = "full",
    ),
    validation_data = (
        id = "1",
        operating_points = PSIDS.SurrogateOperatingPoint[PSIDS.GenerationLoadScale(
            generation_scale = 1.0,
            load_scale = 1.0,
        ),],
        perturbations = repeat(
            [[PSIDS.RandomLoadChange(time = 1.0, load_multiplier_range = (0.0, 2.0))]],
            20,
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
            seed = 100,
        ),
    ),
    model_params = MultiDeviceParams(name = "source_1"),
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
generate_validation_data(p)
validation_dataset = Serialization.deserialize(p.validation_data_path)
display(visualize_dataset(validation_dataset))

sys_validation = System(p.surrogate_system_path)
sys_validation_aux = deepcopy(sys_validation)
PSN.add_surrogate_psid!(sys_validation, p.model_params, validation_dataset)
data_collection_location = Serialization.deserialize(p.data_collection_location_path)[2]

compare_plots = visualize_loss(
    sys_validation,
    sys_validation_aux,
    p_start,
    validation_dataset,
    p.validation_data,
    data_collection_location,
    p.model_params,
)
for p in compare_plots
    display(p)
end

df_loss = evaluate_loss(
    sys_validation,
    sys_validation_aux,
    p_start,
    validation_dataset,
    p.validation_data,
    data_collection_location,
    p.model_params,
)
p1 = scatter(
    df_loss["max_error_ir"],
    xlabel = "fault number",
    ylabel = "Current p.u.",
    label = "max error - Ir",
)
scatter!(
    p1,
    df_loss["max_error_ii"],
    xlabel = "fault number",
    ylabel = "Current p.u.",
    label = "max error - Ii",
)
p2 = scatter(
    df_loss["mae_ir"],
    xlabel = "fault number",
    ylabel = "Current p.u.",
    label = "mae - Ir",
)
scatter!(
    p2,
    df_loss["mae_ii"],
    xlabel = "fault number",
    ylabel = "Current p.u.",
    label = "mae - Ii",
)
display(plot(p1, p2))
