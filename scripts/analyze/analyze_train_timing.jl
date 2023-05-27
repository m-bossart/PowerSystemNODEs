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
#train_folder = joinpath("transfers", "exp_04_25_23_timing")
#train_folder = joinpath("transfers", "exp_05_02_23_timing_nparams")
train_folder = joinpath("transfers", "exp_05_02_23_data_random")
train_folder = joinpath("transfers", "exp_05_11_23_data_random")
train_folder = joinpath("transfers", "exp_05_22_23_data_grid")
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
traces_5 = GenericTrace{Dict{Symbol, Any}}[]
traces_6 =  GenericTrace{Dict{Symbol, Any}}[]
for file in train_files_with_output[2:end]
    @error file
    rebase_path!(file, train_folder)
    params = TrainParams(file)
    path_to_output = joinpath(params.output_data_path, params.train_id)
    output_dict =
        JSON3.read(read(joinpath(path_to_output, "high_level_outputs")), Dict{String, Any})
    n_recorded_iterations = length(output_dict["recorded_iterations"])
    df_loss =
        PowerSimulationNODE.read_arrow_file_to_dataframe(joinpath(path_to_output, "loss"))
    times = df_loss[!, :iteration_time_seconds]
    loss = df_loss[!, :Loss]
    loss_init = df_loss[!, :Loss_initialization]
    loss_dynamic = df_loss[!, :Loss_dynamic]
    df_predictions = PowerSimulationNODE.read_arrow_file_to_dataframe(
        joinpath(path_to_output, "predictions"),
    )
    first_sol = df_predictions[1, "surrogate_solution"]
    deq_iterations = [x.deq_iterations for x in df_predictions[!, "surrogate_solution"]]
    dyn_iterations_accept = []
    dyn_iterations_reject = []
    @show df_predictions[!, "surrogate_solution"][1].destats
    for x in df_predictions[!, "surrogate_solution"]
        if x.destats === nothing
            push!(dyn_iterations_accept, 0.0)
            push!(dyn_iterations_reject, 0.0)
        else
            push!(dyn_iterations_accept, x.destats.naccept)
            push!(dyn_iterations_reject, x.destats.nreject)
        end
    end
    avg_deq_iterations = Statistics.mean(deq_iterations)
    @show dyn_iterations_accept
    @show dyn_iterations_reject
    avg_dyn_iterations_accept = Statistics.mean(filter(x -> x != 0, dyn_iterations_accept))
    #avg_dyn_iterations_reject = Statistics.mean(filter(x -> x != 0, dyn_iterations_reject))
    dynamic_solver = params.optimizer[1].dynamic_solver
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
    push!(
        traces_5,
        PlotlyJS.scatter(
            x = [output_dict["n_params_surrogate"]],
            y = [avg_dyn_iterations_accept],
            name = "$(params.train_id)",
        ),
    )
    push!(
        traces_6,
        PlotlyJS.scatter(
            x = 1:length(times),
            y = loss,
            name = "$(params.train_id) - total",
        ),
    )
    push!(
        traces_6,
        PlotlyJS.scatter(
            x = 1:length(times),
            y = loss_init,
            name = "$(params.train_id) - initialization",
        ),
    )
    push!(
        traces_6,
        PlotlyJS.scatter(
            x = 1:length(times),
            y = loss_dynamic,
            name = "$(params.train_id) - dynamic",
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
p5 = PlotlyJS.plot(
    traces_5,
    Layout(
        xaxis = attr(title = "# of parameters"),
        yaxis = attr(title = "avg accepted dyn iterations"),
    ),
)
p6 = PlotlyJS.plot(
    traces_6,
    Layout(
        xaxis = attr(title = "iteration"),
        yaxis = attr(title = "loss"),
    ),
)

display(p1)
display(p2)
display(p3)
display(p4)
display(p5)
display(p6)

