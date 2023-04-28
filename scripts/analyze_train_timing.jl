using PowerSystems
using PowerSimulationNODE
using PlotlyJS
using JSON3
using Logging
using Serialization
using Statistics

#train_folder = joinpath("transfers", "exp_04_16_23_timing_nparams")
#train_folder = joinpath("transfers", "exp_04_16_23_timing_collocation")
#train_folder = joinpath("transfers", "exp_04_16_23_timing_tstops")
#train_folder = joinpath("transfers", "exp_04_19_23_timing")
#train_folder = joinpath("transfers", "exp_04_21_23_timing")
#train_folder = joinpath("transfers", "exp_04_24_23_timing")
train_folder = joinpath("transfers", "exp_04_25_23_timing")
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
    train_files,
)

traces_1 = GenericTrace{Dict{Symbol, Any}}[]
traces_2 = GenericTrace{Dict{Symbol, Any}}[]
traces_3 = GenericTrace{Dict{Symbol, Any}}[]
traces_4 = GenericTrace{Dict{Symbol, Any}}[]

for file in train_files_with_output
    rebase_path!(file, train_folder)
    params = TrainParams(file)
    path_to_output = joinpath(params.output_data_path, params.train_id)
    output_dict =
        JSON3.read(read(joinpath(path_to_output, "high_level_outputs")), Dict{String, Any})
    n_recorded_iterations = length(output_dict["recorded_iterations"])
    df_loss =
        PowerSimulationNODE.read_arrow_file_to_dataframe(joinpath(path_to_output, "loss"))
    times = df_loss[!, :iteration_time_seconds]
    df_predictions = PowerSimulationNODE.read_arrow_file_to_dataframe(
        joinpath(path_to_output, "predictions"),
    )
    first_sol = df_predictions[1, "surrogate_solution"]
    @error fieldnames(typeof(first_sol))
    deq_iterations = [x.deq_iterations for x in df_predictions[!, "surrogate_solution"]]
    avg_deq_iterations = Statistics.mean(deq_iterations)
    dynamic_solver = params.dynamic_solver
    push!(
        traces_1,
        PlotlyJS.scatter(
            x = [length(first_sol.t_series)],
            y = [Statistics.mean(times[2:end] .- times[1:(end - 1)])],
            name = "$(params.train_id)",
        ),
    )
    push!(
        traces_2,
        PlotlyJS.scatter(
            x = 1:(length(times) - 1),
            y = times[2:end] .- times[1:(end - 1)],
            name = "$(params.train_id)",
        ),
    )
    push!(
        traces_3,
        PlotlyJS.scatter(
            x = [output_dict["n_params_surrogate"]],
            y = [Statistics.mean(times[2:end] .- times[1:(end - 1)])],
            name = "$(params.train_id)",
        ),
    )

    push!(
        traces_4,
        PlotlyJS.scatter(
            x = [output_dict["n_params_surrogate"]],
            y = [avg_deq_iterations],
            name = "$(params.train_id)",
        ),
    )

    #Cleanup to be able to delete the arrow files: https://github.com/apache/arrow-julia/issues/61
    #GC.gc() 
end
p1 = PlotlyJS.plot(
    traces_1,
    Layout(
        xaxis = attr(title = "# of collocation points"),
        yaxis = attr(title = "avg iteration time (s)"),
    ),
)
p2 = PlotlyJS.plot(
    traces_2,
    Layout(
        xaxis = attr(title = "iteration number"),
        yaxis = attr(title = "iteration time (s)"),
    ),
)
p3 = PlotlyJS.plot(
    traces_3,
    Layout(
        xaxis = attr(title = "# of parameters"),
        yaxis = attr(title = "avg iteration time (s)"),
    ),
)
p4 = PlotlyJS.plot(
    traces_4,
    Layout(
        xaxis = attr(title = "# of parameters"),
        yaxis = attr(title = "avg deq iterations"),
    ),
)

display(p1)
display(p2)
display(p3)
display(p4)
