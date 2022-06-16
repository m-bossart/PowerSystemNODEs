#NEW SCRIPT!
#is called in hpc with two arguments: TrainParam file to use for generating data and folder name to save the data in (eg dataset1)
#allows the split of two bash file calls for maximum efficiency. 

using PowerSimulationNODE
using Serialization
using Logging
include("../system_data/dynamic_components_data.jl")

sample_train_parameters = "input_data/sample_parameters.json"
p = TrainParams()
PowerSimulationNODE.serialize(p, "input_data/sample_parameters.json")
train_params_file = isempty(ARGS) ? sample_train_parameters : ARGS[1]
data_folder_name =  isempty(ARGS) ? "dataset_1" : ARGS[2]
train_params = TrainParams(train_params_file)

logger = configure_logging(
    console_level = PowerSimulationNODE.NODE_CONSOLE_LEVEL,
    file_level = PowerSimulationNODE.NODE_FILE_LEVEL,
    filename = "log_generatedata.log",
)

try
    with_logger(logger) do
        train_data_path = joinpath(train_params.input_data_path, "data")
        train_system_path = joinpath(train_params.input_data_path, "system.json")
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
