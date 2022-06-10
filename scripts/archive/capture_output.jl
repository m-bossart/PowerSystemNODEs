
#SCRIPT TO DEBUG ISSUES ON HPC WITH OUT OF MEMORY DUE TO DATA CAPTURE 

using Revise
using PowerSimulationNODE
import DataFrames
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
#Don't get log file written if the process is killed externally... 

function capture_fake_outputs(train_params)
    output_dict = Dict{String, Any}(
        "loss" => DataFrames.DataFrame(
            PVS_name = Vector{String}[],
            RangeCount = Int[],
            Loss = Float64[],
        ),
        "parameters" => DataFrames.DataFrame(Parameters = Vector{Any}[]),
        "predictions" => DataFrames.DataFrame(
            t_prediction = Vector{Any}[],
            prediction = Vector{Any}[],
            observation = Vector{Any}[],
        ),
        "total_time" => [],
        "total_iterations" => 0,
        "recorded_iterations" => [],
        "final_loss" => [],
        "timing_stats_compile" => [],
        "timing_stats" => [],
        "n_params_nn" => 0,
        "train_id" => "",
    )
    output_dict["parameters"] = DataFrames.DataFrame(rand(10000, 2000), :auto)
    output_dict["predictions"] = DataFrames.DataFrame(rand(10000, 2000), :auto)
    output_dict["loss"] = DataFrames.DataFrame(rand(10000, 2000), :auto)

    @info sizeof(output_dict["loss"])
    @time PowerSimulationNODE._capture_output(
        output_dict,
        train_params.output_data_path,
        train_params.train_id,
    )
    return true
end

try
    with_logger(logger) do
        status = capture_fake_outputs(train_params)
    end
    #Add catch ? 
finally
    close(logger)
end
