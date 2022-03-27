using PowerSimulationNODE
using Plots
using Logging
include("../system_data/dynamic_components_data.jl")

configure_logging(console_level = Logging.Info)
#configure_logging(;filename = "train_node.log")

visualize_level = isempty(ARGS) ? 3 : parse(Int64, ARGS[1])

train_files = filter(
    x -> occursin("train_", x),
    readdir(joinpath(pwd(), PowerSimulationNODE.INPUT_FOLDER_NAME), join = true),
)

for file in train_files
    visualize_training(file, visualize_level = visualize_level)
    #animate_training(file, skip_frames = 50, fps = 10)
end
