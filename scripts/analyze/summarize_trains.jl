using PowerSimulationNODE
using Plots

#a = generate_summary(TrainParams().output_data_path)
a = generate_summary(joinpath(pwd(), "exp_1", "output_data"))
#plotlyjs()

p = visualize_summary(a)
print_high_level_output_overview(a, pwd())
png(p, joinpath(TrainParams().base_path, "train_summary"))

ids = []
times = []
for (k, v) in a
    push!(ids, parse(Int64, k))
    push!(times, v["timing_stats"][1]["time"])
end
