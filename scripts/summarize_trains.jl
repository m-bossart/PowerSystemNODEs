using PowerSimulationNODE
using Plots

a = generate_summary(TrainParams().output_data_path)
p = visualize_summary(a)
png(p, joinpath(TrainParams().base_path, "train_summary"))
