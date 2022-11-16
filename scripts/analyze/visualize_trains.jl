using PowerSystems
using PowerSimulationNODE
using Plots
using Logging
#include("../system_data/dynamic_components_data.jl")
train_folder = joinpath("transfers", "exp_11_14_22")
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
##

for file in train_files_with_output
    visualize_training(file, skip = 1, new_base_path = train_folder)       #TODO - re-base path should be separate
    #animate_training(file, skip = 100)
end
