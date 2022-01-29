using PowerSimulationNODE
using Plots

p = visualize_summary(NODETrainParams().output_data_path)
png(p, joinpath(NODETrainParams().base_path, "train_summary"))
