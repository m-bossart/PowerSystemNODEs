function _copy_full_system_to_train_directory(
    SCRATCH_PATH,
    project_folder,
    train_folder,
    system_name,
)
    mkpath(
        joinpath(
            SCRATCH_PATH,
            project_folder,
            train_folder,
            PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
        ),
    )
    cp(
        joinpath(SCRATCH_PATH, project_folder, "systems", string(system_name, ".json")),
        joinpath(
            SCRATCH_PATH,
            project_folder,
            train_folder,
            PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
            string(system_name, ".json"),
        ),
        force = true,
    )
    cp(
        joinpath(
            SCRATCH_PATH,
            project_folder,
            "systems",
            string(system_name, "_validation_descriptors.json"),
        ),
        joinpath(
            SCRATCH_PATH,
            project_folder,
            train_folder,
            PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
            string(system_name, "_validation_descriptors.json"),
        ),
        force = true,
    )
    return
end

function _regenerate_datasets(exp_folder)
    p = TrainParams(joinpath(exp_folder, "input_data", "train_001.json"))
    generate_train_data(p)
    generate_validation_data(p)
    generate_test_data(p)
end

function _get_parameters_from_prior_training(starting_file)
    p_starting = TrainParams(starting_file)
    path_to_output = joinpath(p_starting.output_data_path, p_starting.train_id)
    output_dict =
        JSON3.read(read(joinpath(path_to_output, "high_level_outputs")), Dict{String, Any})
    df_predictions = PowerSimulationNODE.read_arrow_file_to_dataframe(
        joinpath(path_to_output, "predictions"),
    )
    chosen_iteration_index =
        indexin(output_dict["chosen_iteration"], output_dict["recorded_iterations"])[1]
    return df_predictions[chosen_iteration_index, "parameters"][1]
end
function _serialize_starting_parameters(θ)
    mkpath(
        joinpath(
            SCRATCH_PATH,
            project_folder,
            train_folder,
            PowerSimulationNODE.INPUT_FOLDER_NAME,
        ),
    )
end

function determine_p_start(sys, surrogate_buses)
    p_powers = zeros(6)
    p_gfl = zeros(length(PowerSimulationNODE.default_params(PSIDS.GFLParams())))
    p_gfm = zeros(length(PowerSimulationNODE.default_params(PSIDS.GFMParams())))
    p_load = zeros(length(PowerSimulationNODE.default_params(PSIDS.ZIPParams())))

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
            println(
                "name:    ",
                get_name(static_injector),
                "       base power: ",
                get_base_power(static_injector),
            )
            n_gfl += 1
            base_power = get_base_power(static_injector)
            @assert get_dynamic_injector(static_injector) !== nothing
            @assert get_base_power(get_dynamic_injector(static_injector)) ==
                    get_base_power(static_injector)
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:base_power_gfl]] += base_power
            converter = get_converter(i)
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:rated_voltage_gfl]] +=
                get_rated_voltage(converter) * base_power
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:rated_current_gfl]] +=
                get_rated_current(converter) * base_power

            active_power = PowerSystems.get_active_power_control(get_outer_control(i))
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:Kp_p]] +=
                get_Kp_p(active_power) * base_power
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:Ki_p]] +=
                get_Ki_p(active_power) * base_power
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:ωz_gfl]] +=
                get_ωz(active_power) * base_power

            reactive_power = PowerSystems.get_reactive_power_control(get_outer_control(i))
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:Kp_q]] +=
                get_Kp_q(reactive_power) * base_power
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:Ki_q]] +=
                get_Ki_q(reactive_power) * base_power
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:ωf_gfl]] +=
                get_ωf(reactive_power) * base_power

            inner_control = get_inner_control(i)
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:kpc_gfl]] +=
                get_kpc(inner_control) * base_power
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:kic_gfl]] +=
                get_kic(inner_control) * base_power
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:kffv_gfl]] +=
                get_kffv(inner_control) * base_power

            dc_source = get_dc_source(i)
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:voltage_gfl]] +=
                get_voltage(dc_source) * base_power

            freq_estimator = get_freq_estimator(i)
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:ω_lp]] +=
                get_ω_lp(freq_estimator) * base_power
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:kp_pll]] +=
                get_kp_pll(freq_estimator) * base_power
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:ki_pll]] +=
                get_ki_pll(freq_estimator) * base_power

            filter = get_filter(i)
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:lf_gfl]] +=
                get_lf(filter) * base_power
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:rf_gfl]] +=
                get_rf(filter) * base_power
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:cf_gfl]] +=
                get_cf(filter) * base_power
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:lg_gfl]] +=
                get_lg(filter) * base_power
            p_gfl[PowerSimulationNODE.gfl_indices[:params][:rg_gfl]] +=
                get_rg(filter) * base_power
        end
    end
    p_gfl = vcat(p_gfl[1], p_gfl[2:end] ./ p_gfl[1])  #first parameter is base power, every other parameter is base power weighted average.

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
            println(
                "name:    ",
                get_name(static_injector),
                "       base power: ",
                get_base_power(static_injector),
            )
            n_gfm += 1
            base_power = get_base_power(static_injector)
            @assert get_dynamic_injector(static_injector) !== nothing
            @assert get_base_power(get_dynamic_injector(static_injector)) ==
                    get_base_power(static_injector)
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:base_power_gfm]] += base_power
            converter = get_converter(i)
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:rated_voltage_gfm]] +=
                get_rated_voltage(converter) * base_power
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:rated_current_gfm]] +=
                get_rated_current(converter) * base_power

            active_power = PowerSystems.get_active_power_control(get_outer_control(i))
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:Rp]] +=
                get_Rp(active_power) * base_power
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:ωz_gfm]] +=
                get_ωz(active_power) * base_power

            reactive_power = PowerSystems.get_reactive_power_control(get_outer_control(i))
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:kq]] +=
                get_kq(reactive_power) * base_power
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:ωf_gfm]] +=
                get_ωf(reactive_power) * base_power

            inner_control = get_inner_control(i)
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:kpv]] +=
                get_kpv(inner_control) * base_power
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:kiv]] +=
                get_kiv(inner_control) * base_power
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:kffv_gfm]] +=
                get_kffv(inner_control) * base_power
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:rv]] +=
                get_rv(inner_control) * base_power
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:lv]] +=
                get_lv(inner_control) * base_power
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:kpc_gfm]] +=
                get_kpc(inner_control) * base_power
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:kic_gfm]] +=
                get_kpc(inner_control) * base_power
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:kffi]] +=
                get_kffi(inner_control) * base_power
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:ωad]] +=
                get_ωad(inner_control) * base_power
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:kad]] +=
                get_kad(inner_control) * base_power

            dc_source = get_dc_source(i)
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:voltage_gfm]] +=
                get_voltage(dc_source) * base_power

            freq_estimator = get_freq_estimator(i)

            filter = get_filter(i)
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:lf_gfm]] +=
                get_lf(filter) * base_power
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:rf_gfm]] +=
                get_rf(filter) * base_power
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:cf_gfm]] +=
                get_cf(filter) * base_power
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:lg_gfm]] +=
                get_lg(filter) * base_power
            p_gfm[PowerSimulationNODE.gfm_indices[:params][:rg_gfm]] +=
                get_rg(filter) * base_power
        end
    end
    p_gfm = vcat(p_gfm[1], p_gfm[2:end] ./ p_gfm[1])  #first parameter is base power, every other parameter is base power weighted average.

    n_load = 0
    for i in get_components(StandardLoad, sys)
        if get_number(get_bus(i)) in surrogate_buses
            println("name:    ", get_name(i), "       base power: ", get_base_power(i))
            n_load += 1
            base_power = get_base_power(i)
            p_load[PowerSimulationNODE.zip_indices[:params][:base_power_zip]] += base_power
            p_load[PowerSimulationNODE.zip_indices[:params][:max_active_power_Z]] +=
                get_max_impedance_active_power(i) .* base_power
            p_load[PowerSimulationNODE.zip_indices[:params][:max_reactive_power_Z]] +=
                get_max_impedance_reactive_power(i) .* base_power
            p_load[PowerSimulationNODE.zip_indices[:params][:max_active_power_I]] +=
                get_max_current_active_power(i) .* base_power
            p_load[PowerSimulationNODE.zip_indices[:params][:max_reactive_power_I]] +=
                get_max_current_reactive_power(i) .* base_power
            p_load[PowerSimulationNODE.zip_indices[:params][:max_active_power_P]] +=
                get_max_constant_active_power(i) .* base_power
            p_load[PowerSimulationNODE.zip_indices[:params][:max_reactive_power_P]] +=
                get_max_constant_reactive_power(i) .* base_power
        end
    end
    p_load = vcat(p_load[1], p_load[2:end] ./ p_load[1])  #first parameter is base power, every other parameter is base power weighted average.

    p_powers =
        [p_load[1], p_load[1], p_gfl[1], p_gfl[1], p_gfm[1], p_gfm[1]] ./
        (p_load[1] + p_gfl[1] + p_gfm[1]) #hardcoded order... #should this be proportion of base power? 

    p_start = vcat(p_powers, p_load, p_gfl, p_gfm)
    return p_start
end
