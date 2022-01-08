function visualize_training(params::NODETrainParams; visualize_level = 1)
    @debug dump(params)
    path_to_input = joinpath(params.base_path, params.input_data_path)
    path_to_output = joinpath(params.base_path, params.output_data_path, params.train_id)

    output_dict =
        JSON3.read(read(joinpath(path_to_output, "high_level_outputs")), Dict{String, Any})
    println("--------------------------------")
    println("TRAIN ID: ", params.train_id)
    println("TOTAL TIME: ", output_dict["total_time"])
    println("TOTAL ITERATIONS: ", output_dict["total_iterations"])
    println("FINAL LOSS: ", output_dict["final_loss"])
    println("--------------------------------")

    if params.output_mode == 2
        plots = visualize_2(params, path_to_output, path_to_input)
        for p in plots
            png(p, path_to_output)
        end
        return
    elseif params.output_mode == 3
        plots = visualize_3(params, path_to_output, path_to_input, visualize_level)

        for (i, p) in enumerate(plots)
            png(p, joinpath(path_to_output, string("plot_", i)))
        end
        return
    end
end

function visualize_2(params, path_to_output, path_to_input)
    df_loss = DataFrame(Arrow.Table(joinpath(path_to_output, "loss")))
    p1 = plot(df_loss.Loss, title = "Loss")
    p2 = plot(df_loss.RangeCount, title = "Range Count")
    return plot(p1, p2, layout = (2, 1))
end

function visualize_3(params, path_to_output, path_to_input, visualize_level)
    df_loss = DataFrame(Arrow.Table(joinpath(path_to_output, "loss")))
    list_plots = []
    p1 = plot(df_loss.Loss, title = "Loss")
    p2 = plot(df_loss.RangeCount, title = "Range Count")
    p = plot(p1, p2, layout = (2, 1))
    push!(list_plots, p)

    output_dict =
        JSON3.read(read(joinpath(path_to_output, "high_level_outputs")), Dict{String, Any})
    PVS_name = df_loss.PVS_name[:]
    if visualize_level == 1
        transition_indices = [length(PVS_name)]   #LEVEL 1: plot end result only (the result used to calculate final loss)
    elseif visualize_level == 2
        transition_indices = find_transition_indices(PVS_name)  #LEVEL 2: plot when moving to new fault(s)
    elseif visualize_level == 3
        transition_indices = find_transition_indices(df_loss.RangeCount)  #LEVEL 3: plot when moving to new data range
    elseif visualize_level == 4
        transition_indices = collect(1:length(PVS_name)) #LEVEL 4: plot every iteration moving to new data range
    else
        @warn "Invalid value for parameter visualize_level"
    end
    df_predictions = DataFrame(Arrow.Table(joinpath(path_to_output, "predictions")))
    TrainInputs =
        JSON3.read(read(joinpath(params.input_data_path, "data.json")), NODETrainInputs)
    tsteps = TrainInputs.tsteps
    fault_data = TrainInputs.fault_data

    for i in transition_indices
        ir_pred = df_predictions[i, "ir_prediction"]
        ii_pred = df_predictions[i, "ii_prediction"]
        t_pred = df_predictions[i, "t_prediction"]
        i_true = concatonate_i_true(fault_data, df_loss[i, :PVS_name], :)
        ir_true = i_true[1, :]
        ii_true = i_true[2, :]
        t_all = vec(Float64.(concatonate_t(tsteps, df_loss[i, :PVS_name], :))) #TODO fix this syntax 
        p3 = scatter(t_pred, ir_pred, ms = 2, msw = 0, label = "prediction")
        scatter!(p3, t_all, ir_true, ms = 2, msw = 0, label = "truth")
        p4 = scatter(t_pred, ii_pred, ms = 2, msw = 0, label = "prediction")
        scatter!(p4, t_all, ii_true, ms = 2, msw = 0, label = "truth")
        p = plot(
            p3,
            p4,
            title = string(df_loss[i, :PVS_name], " loss: ", output_dict["final_loss"]),
            layout = (2, 1),
        )
        push!(list_plots, p)
    end
    return list_plots
end

function find_transition_indices(list)
    transition_indices = Int[]

    for i in 1:(length(list) - 1)
        if list[i] != list[i + 1]
            push!(transition_indices, i)
        end
    end
    push!(transition_indices, length(list))
    return transition_indices
end

function visualize_summary(output_data_path)
    output_directories = readdir(output_data_path)
    high_level_outputs_dict = Dict{String, Dict{String, Any}}()
    for dir in output_directories
        output_dict = JSON3.read(
            read(joinpath(output_data_path, dir, "high_level_outputs")),
            Dict{String, Any},
        )
        high_level_outputs_dict[output_dict["train_id"]] = output_dict
    end
    p = scatter()
    for (key, value) in high_level_outputs_dict
        scatter!(
            p,
            (value["total_time"], value["final_loss"]),
            label = value["train_id"],
            xlabel = "total time (s)",
            ylabel = "final loss",
            markersize = 3,
            markerstrokewidth = 0,
        )
        annotate!(
            value["total_time"],
            value["final_loss"],
            text(value["train_id"], :red, 8),
        )
    end
    return p
end

#= using Plots, Random
vals = rand(10,2)
p = scatter(vals[:,1], vals[:,2],xlim=[0,1.1])
some_labels=randstring.(fill(5,10))
annotate!.(vals[:,1].+0.01, vals[:,2], text.(some_labels, :red, :left,11))
p =#
