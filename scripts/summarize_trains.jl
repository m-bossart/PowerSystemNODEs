#Script to run after a collection of training instances to provide high level results
include("../src/PowerSystemNODEs.jl")
include("../src/constants.jl")
include("../src/HPCTrain.jl")

p = visualize_summary(NODETrainParams().output_data_path)
png(p, joinpath(NODETrainParams().output_data_path, "train_summary"))
