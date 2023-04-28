using PowerSimulationNODE
using Plots

train_folder = joinpath("transfers", "exp_04_21_23_timing")
plotlyjs()

a = generate_summary(joinpath(pwd(), train_folder, "output_data"))
print_train_parameter_overview(
    joinpath(pwd(), train_folder, PowerSimulationNODE.INPUT_FOLDER_NAME),
)
p = visualize_summary(a)
