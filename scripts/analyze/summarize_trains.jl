using PowerSimulationNODE
using Plots

train_folder = joinpath("transfers", "exp_03_25_23_data_grid")
plotlyjs()  #interactive plot for zooming

a = generate_summary(joinpath(pwd(), train_folder, "output_data"))
print_train_parameter_overview(
    joinpath(pwd(), train_folder, PowerSimulationNODE.INPUT_FOLDER_NAME),
)
p = visualize_summary(a)
##
print_high_level_output_overview(
    a,
    joinpath(pwd(), train_folder, PowerSimulationNODE.INPUT_FOLDER_NAME),
)
png(p, joinpath(TrainParams().base_path, "train_summary"))

ids = []
times = []
for (k, v) in a
    push!(ids, parse(Int64, k))
    push!(times, v["timing_stats"][1]["time"])
end
