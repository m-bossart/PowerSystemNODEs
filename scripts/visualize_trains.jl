include("../src/PowerSystemNODEs.jl")
include("../src/constants.jl")
include("../system_data/dynamic_components_data.jl")
configure_logging(console_level = Logging.Info)
#configure_logging(;filename = "train_node.log")

visualize_level = isempty(ARGS) ? 1 : parse(Int64, ARGS[1])

train_files = filter(
    x -> occursin("train_", x),
    readdir(joinpath(pwd(), INPUT_FOLDER_NAME), join = true),
)
for file in train_files
    p = NODETrainParams(file)
    visualize_training(p, visualize_level = visualize_level)
end
