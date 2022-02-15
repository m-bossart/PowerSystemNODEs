using PowerSimulationNODE
using Serialization
using Logging

include("../system_data/dynamic_components_data.jl")
logger = configure_logging(
    console_level = PowerSimulationNODE.NODE_CONSOLE_LEVEL,
    file_level = PowerSimulationNODE.NODE_FILE_LEVEL,
    filename = "log_generatedata.txt",
)

try
    with_logger(logger) do
        train_data_path = joinpath(PowerSimulationNODE.INPUT_FOLDER_NAME, "data")
        train_system_path = joinpath(PowerSimulationNODE.INPUT_FOLDER_NAME, "system.json")
        full_system_path =
            joinpath(PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME, "full_system.json")
        SURROGATE_BUS = 16

        @warn "Rebuilding full system"
        include("build_full_system.jl")

        @warn "Rebuilding input data files"
        sys_full = node_load_system(full_system_path)
        pvs_data = fault_data_generator("scripts/config.yml", full_system_path)
        sys_pvs = build_pvs(pvs_data)
        label_area!(sys_full, [SURROGATE_BUS], "surrogate")
        @assert check_single_connecting_line_condition(sys_full)
        sys_surr = remove_area(sys_full, "1")
        sys_train = build_train_system(sys_surr, sys_pvs, "surrogate")
        to_json(
            sys_train,
            joinpath(PowerSimulationNODE.INPUT_FOLDER_NAME, "system.json"),
            force = true,
        )
        d = generate_train_data(
            sys_train,
            NODETrainDataParams(ode_model = "none"),
            SURROGATE_BUS,
            inv_case78("aa"),
        )

        Serialization.serialize(train_data_path, d)
    end
finally
    close(logger)
end
######### POST TRAIN GENERATE PREDICTION DATA ########
#= sys_rest = remove_area(sys_full, "surrogate")
sys_reduced = build_reduced_system(sys_rest, NODE, "1")
prediction_data = fault_data_generator(sys_reduced) 
final_loss = loss(ground_truth_data, prediction_data)=#
