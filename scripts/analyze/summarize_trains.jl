using PowerSimulationNODE
using Plots

#a = generate_summary(TrainParams().output_data_path)
a = generate_summary(joinpath(pwd(), "transfers", "exp_11_03_22", "output_data"))
#plotlyjs()
print_train_parameter_overview(joinpath(pwd(), "transfers", "exp_11_03_22", "input_data"))
p = visualize_summary(a)
print_high_level_output_overview(
    a,
    joinpath(pwd(), "transfers", "exp_11_03_22", "input_data"),
)
png(p, joinpath(TrainParams().base_path, "train_summary"))

ids = []
times = []
for (k, v) in a
    push!(ids, parse(Int64, k))
    push!(times, v["timing_stats"][1]["time"])
end
