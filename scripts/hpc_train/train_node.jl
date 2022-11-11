using Revise
using PowerSystems
using PowerSimulationNODE
using Logging

train_params_file = ARGS[1]
train_params = TrainParams(train_params_file)

logger = configure_logging(
    console_level = PowerSimulationNODE.NODE_CONSOLE_LEVEL,
    file_level = PowerSimulationNODE.NODE_FILE_LEVEL,
    filename = joinpath(
        train_params.base_path,
        string("log_", train_params.train_id, ".log"),
    ),
)
try
    with_logger(logger) do
        @info train_params.train_id
        status = train(train_params)
    end
catch
    close(logger)
end
