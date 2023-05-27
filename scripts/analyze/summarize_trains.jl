using PowerSimulationNODE
using Plots

train_folder = joinpath("transfers", "exp_05_22_23_data_grid")
plotlyjs()

a = generate_summary(joinpath(pwd(), train_folder, "output_data"))
print_train_parameter_overview(
    joinpath(pwd(), train_folder, PowerSimulationNODE.INPUT_FOLDER_NAME),
)
p = visualize_summary(a)
##
iterations = []
complete = 0.0
for (k, v) in a
    println(k, "     ", v["total_iterations"])
    push!(iterations, v["total_iterations"])
    if v["total_iterations"] == 15000
        complete += 1.0
    end
end
Plots.scatter(iterations)
