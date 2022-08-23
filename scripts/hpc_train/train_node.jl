using Revise
using PowerSystems
using PowerSimulationNODE
using Logging

train_params_file = ARGS[1]
train_params = TrainParams(train_params_file)

logger = configure_logging(
    console_level = PowerSimulationNODE.NODE_CONSOLE_LEVEL,
    file_level = PowerSimulationNODE.NODE_FILE_LEVEL,
    filename = string("log_", train_params.train_id, ".log"),
)
#TODO - We don't get log file written if the process is killed externally... 
try
    with_logger(logger) do
        status = train(train_params)
    end
catch
    close(logger)
end
