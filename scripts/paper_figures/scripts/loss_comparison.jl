using Plots
using PowerSimulationNODE
using Serialization
using LaTeXStrings

########### INPUT DATA ########### 
exp_folder_data = "transfers/exp_03_06_23_data_random"
train_id_data = "026"

exp_folder_physics = "transfers/exp_03_02_23_physics_random"
train_id_physics = "002"

plot_iterations = 1:2000
y_scale = (0.0001, 1)

################################# 
df_loss_data = PowerSimulationNODE.read_arrow_file_to_dataframe(
    joinpath(exp_folder_data, "output_data", train_id_data, "loss"),
)
df_loss_physics = PowerSimulationNODE.read_arrow_file_to_dataframe(
    joinpath(exp_folder_physics, "output_data", train_id_physics, "loss"),
)
p1 = plot(
    plot_iterations,
    df_loss_data[plot_iterations, :Loss_dynamic],
    yaxis = :log,
    xlabel = "Iteration",
    ylabel = "Loss",
    ylims = y_scale,
    label = L"$L_{\text{dyn}}$",
    labelfontsize = 12,
    tickfontsize = 10,
)
plot!(
    p1,
    plot_iterations,
    df_loss_data[plot_iterations, :Loss_initialization],
    yaxis = :log,
    ylims = y_scale,
    label = L"$L_{\text{dyn}}$",
)
p2 = plot(
    plot_iterations,
    df_loss_physics[plot_iterations, :Loss_dynamic],
    yaxis = :log,
    xlabel = "Iteration",
    ylabel = "Loss",
    ylims = y_scale,
    label = L"$L_{\text{dyn}}$",
    labelfontsize = 12,
    tickfontsize = 10,
)
p = plot(p1, p2, layout = (2, 1), size = (500, 400))
png(p, "figs/23_paperfigs/loss_comparison.png")
