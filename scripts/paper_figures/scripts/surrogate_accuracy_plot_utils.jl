function show_parameter_change(file)
    params = TrainParams(file)
    path_to_output = joinpath(file, "..", "..", "output_data", params.train_id)
    df_predictions = PowerSimulationNODE.read_arrow_file_to_dataframe(
        joinpath(path_to_output, "predictions"),
    )
    θ = df_predictions[:, "parameters"]
    p = Plots.plot()
    percent_change = (θ[1][1] .- θ[end][1]) ./ θ[1][1]
    Plots.scatter!(1:length(θ[end][1]), percent_change)
    display(p)
    return
end

function generate_and_serialize_dataset(file, dataset_type; chosen = 0)
    params = TrainParams(file)
    path_to_output = joinpath(params.output_data_path, params.train_id)
    output_dict =
        JSON3.read(read(joinpath(path_to_output, "high_level_outputs")), Dict{String, Any})

    validation_sys = node_load_system(params.modified_surrogate_system_path)
    validation_sys_aux = node_load_system(params.surrogate_system_path)
    if dataset_type == "test"
        dataset = Serialization.deserialize(params.test_data_path)
    elseif dataset_type == "validation"
        dataset = Serialization.deserialize(params.validation_data_path)
    elseif dataset_type == "train"
        dataset = Serialization.deserialize(params.train_data_path)
    end

    data_collection_location =
        Serialization.deserialize(params.data_collection_location_path)[2]
    df_predictions = PowerSimulationNODE.read_arrow_file_to_dataframe(
        joinpath(path_to_output, "predictions"),
    )
    chosen_iteration_index =
        indexin(output_dict["chosen_iteration"], output_dict["recorded_iterations"])[1]
    if chosen !== 0
        chosen_iteration_index = chosen
    end
    @warn "chosen_iteration_index = $chosen_iteration_index"
    θ = df_predictions[chosen_iteration_index, "parameters"][1]
    @warn length(θ)
    display(validation_sys)
    PowerSimulationNODE.parameterize_surrogate_psid!(validation_sys, θ, params.model_params)
    display(Simulation(MassMatrixModel, validation_sys, pwd(), (0.0, 1.0)))

    if dataset_type == "test"
        surrogate_dataset = generate_surrogate_dataset(
            validation_sys,
            validation_sys_aux,
            θ,
            dataset,
            params.test_data,
            data_collection_location,
            params.model_params,
        )
    elseif dataset_type == "validation"
        surrogate_dataset = generate_surrogate_dataset(
            validation_sys,
            validation_sys_aux,
            θ,
            dataset,
            params.validation_data,
            data_collection_location,
            params.model_params,
        )
    elseif dataset_type == "train"
        surrogate_dataset = generate_surrogate_dataset(
            validation_sys,
            validation_sys_aux,
            θ,
            dataset,
            params.train_data,
            data_collection_location,
            params.model_params,
        )
    end

    Serialization.serialize(
        joinpath(
            params.output_data_path,
            params.train_id,
            string("surrogate_", dataset_type, "_", string(chosen), "_dataset"),
        ),
        surrogate_dataset,
    )
end

#Change this for only the ones where both datasets are stable
function all_current_errors(dataset_1, dataset_2)
    ir_errors = Float64[]
    ii_errors = Float64[]
    for (x, y) in zip(dataset_1, dataset_2)
        if x.stable && y.stable
            ir_errors = vcat(ir_errors, abs.(x.real_current .- y.real_current)')
            ii_errors = vcat(ii_errors, abs.(x.imag_current .- y.imag_current)')
        end
    end
    return vec(ir_errors), vec(ii_errors)
end

function loss_metrics(dataset_surrogate, dataset_groundtruth)
    times = [s.solve_time for s in dataset_surrogate]
    avg_time = sum(times) / length(times)
    println("all times:", times)
    println("average times:   ", avg_time)

    times_ground_truth = [s.solve_time for s in dataset_groundtruth]
    avg_time_ground_truth = sum(times_ground_truth) / length(times_ground_truth)
    println("all times ground truth:   ", times_ground_truth)
    println("average time ground truth:   ", avg_time_ground_truth)

    loss_dict = evaluate_loss(dataset_surrogate, dataset_groundtruth)
    return loss_dict, times, times_ground_truth
end

function _regenerate_datasets(dataset_to_compare, results_to_compare)
    for r in results_to_compare
        if r.generate_data == true
            file = joinpath(r.exp_folder, "input_data", string("train_", r.train_id, ".json"))
            generate_and_serialize_dataset(
                file,
                dataset_to_compare;
                chosen = r.chosen_iteration,
            )
        end
    end
end 

function _plot_historgram_all_errors(dataset_to_compare, results_to_compare)
    traces_ir = GenericTrace{Dict{Symbol, Any}}[]
    traces_ii = GenericTrace{Dict{Symbol, Any}}[]
    for r in results_to_compare
        file = joinpath(r.exp_folder, "input_data", string("train_", r.train_id, ".json"))
        params = TrainParams(file)
        if dataset_to_compare == "train"
            data_ground_truth = Serialization.deserialize(params.train_data_path)
        elseif dataset_to_compare == "validation"
            data_ground_truth = Serialization.deserialize(params.validation_data_path)
        elseif dataset_to_compare == "test"
            data_ground_truth = Serialization.deserialize(params.test_data_path)
        end
        data_surrogate = Serialization.deserialize(
            joinpath(
                params.output_data_path,
                params.train_id,
                string(
                    "surrogate_",
                    dataset_to_compare,
                    "_",
                    string(r.chosen_iteration),
                    "_dataset",
                ),
            ),
        )
        ir_errors, ii_errors = all_current_errors(data_ground_truth, data_surrogate)
        trace_ir = PlotlyJS.histogram(x = ir_errors, opacity = 0.6, name = r.name)
        trace_ii = PlotlyJS.histogram(x = ii_errors, opacity = 0.6, name = r.name)
        push!(traces_ir, trace_ir)
        push!(traces_ii, trace_ii)
    end
    
    p1 = PlotlyJS.plot(
            traces_ir,
            Layout(
                barmode = "overlay",
                xaxis = attr(font_size = 12),
                yaxis = attr(
                    title = "Distribution of errors in Ir (p.u.)",
                    font_size = 12,
                    zeroline = false,
                ),
                legend = attr(
                    x = 1,
                    y = 1.0,
                    font_size = 12,
                    yanchor = "bottom",
                    xanchor = "right",
                    orientation = "h",
                ),
            ),
        )
    
     p2 = PlotlyJS.plot(
            traces_ii,
            Layout(
                barmode = "overlay",
                xaxis = attr(font_size = 12),
                yaxis = attr(
                    title = "Distribution of errors in Ii (p.u.)",
                    font_size = 12,
                    zeroline = false,
                ),
                legend = attr(
                    x = 1,
                    y = 1.0,
                    font_size = 12,
                    yanchor = "bottom",
                    xanchor = "right",
                    orientation = "h",
                ),
            ),
        )
    return p1, p2 
end 

function _plot_box_plot_mean_errors(dataset_to_compare, results_to_compare)
    traces = GenericTrace{Dict{Symbol, Any}}[]
for r in results_to_compare
    file = joinpath(r.exp_folder, "input_data", string("train_", r.train_id, ".json"))
    params = TrainParams(file)
    if dataset_to_compare == "train"
        data_ground_truth = Serialization.deserialize(params.train_data_path)
    elseif dataset_to_compare == "validation"
        data_ground_truth = Serialization.deserialize(params.validation_data_path)
    elseif dataset_to_compare == "test"
        data_ground_truth = Serialization.deserialize(params.test_data_path)
    end
    data_surrogate = Serialization.deserialize(
        joinpath(
            params.output_data_path,
            params.train_id,
            string(
                "surrogate_",
                dataset_to_compare,
                "_",
                string(r.chosen_iteration),
                "_dataset",
            ),
        ),
    )

    loss_dict, times_surrogate, times_ground_truth, =
        loss_metrics(data_surrogate, data_ground_truth)

    ir_mean_error = loss_dict["mae_ir"]
    ii_mean_error = loss_dict["mae_ii"]

    x0_data =
        vcat(repeat(["Iᵣ"], length(ir_mean_error)), repeat(["Iᵢ"], length(ii_mean_error)))

    trace1 = box(;
        y = vcat(ir_mean_error, ii_mean_error),
        x = x0_data,
        boxpoints = "all",
        name = r.name,
        legendgroup = "1",
        #marker_color = :blue,
        marker_size = 3,
        #jitter = 100,
        #pointpos = 0,
    )
    push!(traces, trace1)
end
layout = Layout(;
    font_family = "Times New Roman",
    xaxis = attr(font_size = 12,showline=true, linewidth=1, linecolor="black"),
    yaxis = attr(title = "MAE per fault (p.u. current)",font_size = 12, zeroline = false, showline=true, linewidth=1, linecolor="black"),
    legend = attr(
        x = 1,
        y = 1.0,
        font_size = 12,
        yanchor = "bottom",
        xanchor = "right",
        orientation = "h",
    ),
    #width = 50,
    #height = 40,
    boxmode = "group",
    template = "plotly_white",
    margin = (b =20, t=0, r = 50, l =50)
)

return PlotlyJS.plot(traces, layout)
end 

function  _plot_timing_comparison(dataset_to_compare, results_to_compare)
    traces = GenericTrace{Dict{Symbol, Any}}[]
times_ground_truth = nothing
for r in results_to_compare
    file = joinpath(r.exp_folder, "input_data", string("train_", r.train_id, ".json"))
    params = TrainParams(file)
    if dataset_to_compare == "train"
        data_ground_truth = Serialization.deserialize(params.train_data_path)
    elseif dataset_to_compare == "validation"
        data_ground_truth = Serialization.deserialize(params.validation_data_path)
    elseif dataset_to_compare == "test"
        data_ground_truth = Serialization.deserialize(params.test_data_path)
    end
    data_surrogate = Serialization.deserialize(
        joinpath(
            params.output_data_path,
            params.train_id,
            string(
                "surrogate_",
                dataset_to_compare,
                "_",
                string(r.chosen_iteration),
                "_dataset",
            ),
        ),
    )

    loss_dict, times_surrogate, times_ground_truth, =
        loss_metrics(data_surrogate, data_ground_truth)
    x0_data = repeat([""], length(times_surrogate))
    times_surrogate_stable = filter(x -> x !== 0.0, times_surrogate)
    times_ground_truth_stable = filter(x -> x !== 0.0, times_ground_truth)

    surrogate_stable = [d.stable for d in data_surrogate]
    groundtruth_stable = [d.stable for d in data_ground_truth]

    tp = [x == 1 && y == 1 for (x, y) in zip(surrogate_stable, groundtruth_stable)]
    tn = [x == 0 && y == 0 for (x, y) in zip(surrogate_stable, groundtruth_stable)]
    fp = [x == 1 && y == 0 for (x, y) in zip(surrogate_stable, groundtruth_stable)]
    fn = [x == 0 && y == 1 for (x, y) in zip(surrogate_stable, groundtruth_stable)]
    println("True positives: ", sum(tp))
    println("True negatives: ", sum(tn))
    println("False positives: ", sum(fp))
    println("False negatives: ", sum(fn))

    println(r.exp_folder)
    println(length(times_ground_truth_stable))
    println(length(times_surrogate_stable))
    global times_gt = times_ground_truth_stable
    trace1 = box(;
        y = times_surrogate_stable,
        x = x0_data,
        boxpoints = "all",
        name = r.name,
        legendgroup = "1",
        #marker_color = :blue,
        marker_size = 3,
        #jitter = 100,
        #pointpos = 0,
    )
    push!(traces, trace1)
end
x0_data = repeat([""], length(times_gt))

trace1 = box(;
    y = times_gt,
    x = x0_data,
    boxpoints = "all",
    name = "ground truth",
    legendgroup = "1",
    marker_color = :black,
    marker_size = 3,
    #jitter = 100,
    #pointpos = 0,
)
push!(traces, trace1)

layout = Layout(;
    font_family = "Times New Roman",
    xaxis = attr(font_size = 12, showline=true, linewidth=1, linecolor="black"),
    yaxis = attr(title = "Simulation time per fault (s)", font_size = 12, zeroline = false, showline=true, linewidth=1, linecolor="black"),
    legend = attr(
        x = 1,
        y = 1.0,
        font_size = 12,
        yanchor = "bottom",
        xanchor = "right",
        orientation = "h",
    ),
    boxmode = "group",
    template = "plotly_white",
    margin = (b =20, t=30, r = 50, l =50)
)
return  PlotlyJS.plot(traces, layout)
end 




function _display_comparisons_individual_traces(dataset_to_compare, results_to_compare)
    f1 = joinpath(
        results_to_compare[1].exp_folder,
        "input_data",
        string("train_", results_to_compare[1].train_id, ".json"),
    )
    p1 = TrainParams(f1)

    if dataset_to_compare == "train"
        data_ground_truth = Serialization.deserialize(p1.train_data_path)
    elseif dataset_to_compare == "validation"
        data_ground_truth = Serialization.deserialize(p1.validation_data_path)
    elseif dataset_to_compare == "test"
        data_ground_truth = Serialization.deserialize(p1.test_data_path)
    end

    for (ix, d) in enumerate(data_ground_truth)
        if d.stable == true 
            traces_2 = GenericTrace{Dict{Symbol, Any}}[] 
            for r in results_to_compare 
                file =  joinpath(r.exp_folder, "input_data", string("train_", r.train_id, ".json"))
                params = TrainParams(file)
                dataset_surrogate = Serialization.deserialize(joinpath(params.output_data_path, params.train_id, string("surrogate_", dataset_to_compare, "_", string(r.chosen_iteration), "_dataset")))
                trace_surrogate = PlotlyJS.scatter(;x=dataset_surrogate[ix].tsteps, y=vec(dataset_surrogate[ix].real_current), mode="lines",  name =r.name)
                push!(traces_2, trace_surrogate)
            end 
            trace_ground_truth =  PlotlyJS.scatter(;x=d.tsteps, y=vec(d.real_current), mode="lines", name = "ground truth", line_color = "black")
            #trace_ground_truth =  PlotlyJS.scatter(;x=[1,2], y=[3,4], mode="lines", name = "ground truth")
            push!(traces_2, trace_ground_truth)
            display(PlotlyJS.plot(traces_2))
        end 
    end  
end 