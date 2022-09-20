using PowerSystems
using PowerSimulationNODE
using Serialization
using Logging

include("../system_data/dynamic_components_data.jl")
logger = configure_logging(
    console_level = PowerSimulationNODE.NODE_CONSOLE_LEVEL,
    file_level = PowerSimulationNODE.NODE_FILE_LEVEL,
    filename = "log_generatedata.log",
)

try
    with_logger(logger) do
        train_data_path = joinpath(PowerSimulationNODE.INPUT_FOLDER_NAME, "data")
        train_system_path = joinpath(PowerSimulationNODE.INPUT_FOLDER_NAME, "system.json")
        full_system_path =
            joinpath(PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME, "full_system.json")
        SURROGATE_BUS = 16

        include("build_full_system.jl") #TODO make this more flexible, include multiple systems

        sys_full = node_load_system(full_system_path)
        label_area!(sys_full, [SURROGATE_BUS], "surrogate")

        pvs_coeffs = Dict{Int, Array{NamedTuple}}()
        pvs_coeffs[1] = [(
            internal_voltage_frequencies = [2 * pi / 3],
            internal_voltage_coefficients = [(0.001, 0.0)],
            internal_angle_frequencies = [2 * pi / 3],
            internal_angle_coefficients = [(0.0, 0.0)],
        )]

        pvs_data = generate_pvs_data(sys_full, pvs_coeffs, "surrogate")
        sys_train = create_surrogate_training_system(sys_full, "surrogate", pvs_data)
        d = generate_train_data(
            sys_train,
            GenerateDataParams(tspan = (0.0, 4.0), steps = 400),
        )

        PSY.to_json(sys_train, train_system_path, force = true)

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
