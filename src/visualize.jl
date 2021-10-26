function visualize_training(params::NODETrainParams)
    dump(params)
    df_output = DataFrame(Arrow.Table(joinpath(params.export_file_path, "outputdata")))
    df_input = DataFrame(Arrow.Table(data.faults_results_dir))
    if params.export_mode == 1
        return visualize_1(df_input, df_output)
    elseif params.export_mode == 2
        return visualize_2(df_input, df_output)
    elseif params.export_mode == 3
        return visualize_3(df_input, df_output)
    end
end

function visualize_1(df_input, df_output)
    println("visualize 1 ")
    println("FINAL LOSS: ", df_output[end, 2])
    println("FINAL PARAMETERS: ", df_output[end, 1])
end

function visualize_2(df_input, df_output)
    p1 = plot(df_output[:, 3], title = "Loss")
    p2 = plot(df_output[:, 2], title = "Range count")
    p = plot(p1, p2, layout = (2, 1))
    println("FINAL LOSS: ", df_output[end, 3])
    return p
end

function visualize_3(df_input, df_output)
    p1 = plot(df_output[:, 3], title = "Loss")
    p2 = plot(df_output[:, 2], title = "Range count")
    p3 = plot(df_input.ir_true, title = "ir", label = "true")
    plot!(p3, df_output[end, 5], label = "pred")
    p4 = plot(df_input.ii_true, title = "ii", label = "true")
    plot!(p4, df_output[end, 6], label = "pred")
    p = plot(p1, p3, p2, p4, layout = (2, 2))
    println("FINAL LOSS: ", df_output[end, 3])
    return p

    #Add option for gif? 
    #=anim = Animation()
    for plt in list_plots
        frame(anim, plt)
    end =#
end
