using PowerSimulationNODE
using Serialization
using Logging

train_params_file = split(ARGS[1], ",")[1]
dataset_type =  split(ARGS[1], ",")[2]
p = TrainParams(train_params_file)

logger = configure_logging(
    console_level = PowerSimulationNODE.NODE_CONSOLE_LEVEL,
    file_level = PowerSimulationNODE.NODE_FILE_LEVEL,
    filename = "log_generatedata.log",
)

try
    with_logger(logger) do
        if dataset_type == "train"
            generate_train_data(p)
        elseif dataset_type == "validation"
            generate_validation_data(p)
        elseif dataset_type == "test"
            generate_test_data(p)
        else
            @error "Invalid type of dataset: must be train, validation, or test"
        end
    end
finally
    close(logger)
end
