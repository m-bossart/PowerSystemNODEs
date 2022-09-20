using PowerSimulationNODE
using Plots

a = generate_summary(TrainParams().output_data_path)
#plotlyjs()

p = visualize_summary(a)
print_high_level_output_overview(a, pwd())
png(p, joinpath(TrainParams().base_path, "train_summary"))
