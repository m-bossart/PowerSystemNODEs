using PowerSystems
using PowerSimulationNODE
using Plots
using JSON3
using Logging
#include("../system_data/dynamic_components_data.jl")
train_folder = joinpath("transfers", "exp_01_19_23")

configure_logging(console_level = Logging.Info)
visualize_level = isempty(ARGS) ? 3 : parse(Int64, ARGS[1])

train_files = filter(
    x -> occursin("train_", x) && occursin(".json", x),
    readdir(
        joinpath(pwd(), train_folder, PowerSimulationNODE.INPUT_FOLDER_NAME),
        join = true,
    ),
)

output_folders = readdir(
    joinpath(pwd(), train_folder, PowerSimulationNODE.OUTPUT_FOLDER_NAME),
    join = false,
)

train_files_with_output = filter(
    x ->
        occursin("train_", x) &&
            occursin(".json", x) &&
            TrainParams(x).train_id in output_folders,
    readdir(
        joinpath(pwd(), train_folder, PowerSimulationNODE.INPUT_FOLDER_NAME),
        join = true,
    ),
)

for file in train_files_with_output
    rebase_path!(file, train_folder)
    params = TrainParams(file)
    path_to_output = joinpath(file, "..", "..", "output_data", params.train_id)
    output_dict =
        JSON3.read(read(joinpath(path_to_output, "high_level_outputs")), Dict{String, Any})
    n_recorded_iterations = length(output_dict["recorded_iterations"])
    visualize_training(file, vcat(2:6, (n_recorded_iterations - 4):n_recorded_iterations))
    #animate_training(file, skip = 100)
end
