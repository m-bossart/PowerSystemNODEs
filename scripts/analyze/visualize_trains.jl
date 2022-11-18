using PowerSystems
using PowerSimulationNODE
using Plots
using Logging
#include("../system_data/dynamic_components_data.jl")
train_folder = joinpath("transfers", "exp_11_16_22")

configure_logging(console_level = Logging.Info)
visualize_level = isempty(ARGS) ? 3 : parse(Int64, ARGS[1])

train_files = filter(
    x -> occursin("train_", x) && occursin(".json", x),
    readdir(
        joinpath(pwd(), train_folder, PowerSimulationNODE.INPUT_FOLDER_NAME),
        join = true,
    ),
)

output_folders = readdir(
    joinpath(pwd(), train_folder, PowerSimulationNODE.OUTPUT_FOLDER_NAME),
    join = false,
)

train_files_with_output = filter(
    x ->
        occursin("train_", x) &&
            occursin(".json", x) &&
            TrainParams(x).train_id in output_folders,
    readdir(
        joinpath(pwd(), train_folder, PowerSimulationNODE.INPUT_FOLDER_NAME),
        join = true,
    ),
)
##

for file in train_files_with_outputs
    visualize_training(file, skip = 100, new_base_path = train_folder)       #TODO - re-base path should be separate
    #animate_training(file, skip = 100)
end

L19 = PowerSimulationNODE.read_arrow_file_to_dataframe(
    joinpath(pwd(), "transfers", "exp_11_16_22", "output_data", "019", "loss"),
)
L21 = PowerSimulationNODE.read_arrow_file_to_dataframe(
    joinpath(pwd(), "transfers", "exp_11_16_22", "output_data", "021", "loss"),
)
L23 = PowerSimulationNODE.read_arrow_file_to_dataframe(
    joinpath(pwd(), "transfers", "exp_11_16_22", "output_data", "023", "loss"),
)
L25 = PowerSimulationNODE.read_arrow_file_to_dataframe(
    joinpath(pwd(), "transfers", "exp_11_16_22", "output_data", "025", "loss"),
)
L27 = PowerSimulationNODE.read_arrow_file_to_dataframe(
    joinpath(pwd(), "transfers", "exp_11_16_22", "output_data", "027", "loss"),
)
L29 = PowerSimulationNODE.read_arrow_file_to_dataframe(
    joinpath(pwd(), "transfers", "exp_11_16_22", "output_data", "029", "loss"),
)
L31 = PowerSimulationNODE.read_arrow_file_to_dataframe(
    joinpath(pwd(), "transfers", "exp_11_16_22", "output_data", "031", "loss"),
)
L33 = PowerSimulationNODE.read_arrow_file_to_dataframe(
    joinpath(pwd(), "transfers", "exp_11_16_22", "output_data", "033", "loss"),
)
L35 = PowerSimulationNODE.read_arrow_file_to_dataframe(
    joinpath(pwd(), "transfers", "exp_11_16_22", "output_data", "035", "loss"),
)
p1 = plot(
    L19[!, :Loss_dynamic],
    yscale = :log,
    title = "dataset seed=1",
    label = "no restart",
)
plot!(p1, L25[!, :Loss_dynamic], label = "1 restart")
plot!(p1, L31[!, :Loss_dynamic], label = "2 restart")
png(p1, joinpath(pwd(), "dataset1"))

p2 = plot(
    L21[!, :Loss_dynamic],
    yscale = :log,
    title = "dataset seed=2",
    label = "no restart",
)
plot!(p2, L27[!, :Loss_dynamic], label = "1 restart")
plot!(p2, L33[!, :Loss_dynamic], label = "2 restart")
png(p2, joinpath(pwd(), "dataset2"))

p3 = plot(
    L23[!, :Loss_dynamic],
    yscale = :log,
    title = "dataset seed=3",
    label = "no restart",
)
plot!(p3, L29[!, :Loss_dynamic], label = "1 restart")
plot!(p3, L35[!, :Loss_dynamic], label = "2 restart")
png(p3, joinpath(pwd(), "dataset3"))

p4 = plot(
    L23[1:end, :Loss_dynamic],
    yscale = :log,
    legend = :topleft,
    label = "dynamic loss",
)
plot!(
    p4,
    L23[1:end, :Loss_initialization],
    yscale = :log,
    legend = :topleft,
    label = "initialization loss",
)
png(p4, joinpath(pwd(), "adam_bfgs"))
p4 = plot(
    L23[2900:end, :Loss_dynamic],
    yscale = :log,
    legend = :topleft,
    label = "dynamic loss",
)
plot!(
    p4,
    L23[2900:end, :Loss_initialization],
    yscale = :log,
    legend = :topleft,
    label = "initialization loss",
)
png(p4, joinpath(pwd(), "adam_bfgs_zoom"))
