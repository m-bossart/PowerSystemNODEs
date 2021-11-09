function visualize_training(params::NODETrainParams)
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
        return visualize_2(params, path_to_output, path_to_input)
    elseif params.output_mode == 3
        return visualize_3(params, path_to_output, path_to_input)
    end
end

function visualize_2(params, path_to_output, path_to_input)
    df_loss = DataFrame(Arrow.Table(joinpath(path_to_output, "loss")))
    p1 = plot(df_loss.Loss, title = "Loss")
    p2 = plot(df_loss.RangeCount, title = "Range Count")
    return plot(p1, p2, layout = (2, 1))
end

function visualize_3(params, path_to_output, path_to_input)
    df_loss = DataFrame(Arrow.Table(joinpath(path_to_output, "loss")))
    list_plots = []
    p1 = plot(df_loss.Loss, title = "Loss")
    p2 = plot(df_loss.RangeCount, title = "Range Count")
    p = plot(p1, p2, layout = (2, 1))
    push!(list_plots, p)

    output_dict =
        JSON3.read(read(joinpath(path_to_output, "high_level_outputs")), Dict{String, Any})
    IDs = df_loss.ID[:]
    transition_indices = find_transition_indices(IDs)
    df_predictions = DataFrame(Arrow.Table(joinpath(path_to_output, "predictions")))
    input_dict = JSON3.read(
        read(joinpath(path_to_input, "data.json")),
        Dict{String, Dict{Symbol, Any}},
    )
    for i in transition_indices
        ir_pred = df_predictions[i, "ir_prediction"]
        ii_pred = df_predictions[i, "ii_prediction"]
        ir_true = input_dict[IDs[i]][:ir_ground_truth]
        ii_true = input_dict[IDs[i]][:ii_ground_truth]
        p3 = plot(ir_pred, label = "prediction")
        plot!(p3, ir_true, label = "truth")
        p4 = plot(ii_pred, label = "prediction")
        plot!(p4, ii_true, label = "truth")
        p = plot(
            p3,
            p4,
            title = string(IDs[i], " loss: ", output_dict["final_loss"]),
            layout = (2, 1),
        )
        push!(list_plots, p)
    end
    return list_plots
end

function find_transition_indices(list)
    transition_indices = Int[]

    for i in 1:(length(list) - 1)
        if list[i] !== list[i + 1]
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
        )
    end
    return p
end
