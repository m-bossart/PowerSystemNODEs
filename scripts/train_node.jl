using Revise
using PowerSimulationNODE
using Plots
using Logging
include("../system_data/dynamic_components_data.jl")

sample_train_parameters = "input_data/sample_parameters.json"
p = NODETrainParams()
p.verify_psid_node_off = false
PowerSimulationNODE.serialize(p, "input_data/sample_parameters.json")
train_params_file = isempty(ARGS) ? sample_train_parameters : ARGS[1]
train_params = NODETrainParams(train_params_file)

logger = configure_logging(
    console_level = PowerSimulationNODE.NODE_CONSOLE_LEVEL,
    file_level = PowerSimulationNODE.NODE_FILE_LEVEL,
    filename = string("log_", train_params.train_id, ".log"),
)
try
    with_logger(logger) do
        status = train(train_params)
    end
finally
    close(logger)
end
