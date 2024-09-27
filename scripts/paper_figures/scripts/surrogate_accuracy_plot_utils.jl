using Statistics 
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
    if chosen == 0 
        chosen_iteration_index =
            indexin(output_dict["chosen_iteration"], output_dict["recorded_iterations"])[1]
    else chosen !== 0
        chosen_iteration_index = indexin(chosen, output_dict["recorded_iterations"])[1]
    end
    @warn "chosen_iteration_index = $chosen_iteration_index"
    θ = df_predictions[chosen_iteration_index, "parameters"][1]
    @warn length(θ)
    display(validation_sys)
    PowerSimulationNODE.parameterize_surrogate_psid!(validation_sys, θ, params.model_params)
    display(Simulation(MassMatrixModel, validation_sys, pwd(), (0.0, 1.0)))

    if dataset_type == "test"
        dataset_id = params.test_data.id
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
        dataset_id = params.validation_data.id
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
        dataset_id = params.train_data.id
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
            string("surrogate_", dataset_type, "_", string(dataset_id), "_", string(chosen), "_dataset"),
        ),
        surrogate_dataset,
    )
end

function all_current_errors(dataset_1, dataset_2)
    ir_errors = Float64[]
    ii_errors = Float64[]
    for (x, y) in zip(dataset_1, dataset_2)
        if x.stable && y.stable
            ir_errors = vcat(ir_errors, abs.(get_device_terminal_data(x)[:ir] .- get_device_terminal_data(y)[:ir]))
            ii_errors = vcat(ii_errors, abs.(get_device_terminal_data(x)[:ii] .- get_device_terminal_data(y)[:ii]))
        end
    end
    return vec(ir_errors), vec(ii_errors)
end

function loss_metrics(dataset_surrogate, dataset_groundtruth)
    times = [s.solve_time for s in dataset_surrogate]
    times = filter(x-> x != 0.0, times)
    avg_time = sum(times) / length(times)
    println("average times:   ", avg_time)

    times_ground_truth = [s.solve_time for s in dataset_groundtruth]
    times_ground_truth = filter(x-> x != 0.0, times_ground_truth)
    avg_time_ground_truth = sum(times_ground_truth) / length(times_ground_truth)
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
            dataset_id = params.train_data.id
            data_ground_truth = Serialization.deserialize(params.train_data_path)
        elseif dataset_to_compare == "validation"
            dataset_id = params.validation_data.id
            data_ground_truth = Serialization.deserialize(params.validation_data_path)
        elseif dataset_to_compare == "test"
            dataset_id = params.test_data.id
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
                    string(dataset_id),
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
                    x = -0.5,
                    y = 1.0,
                    font_size = 12,
                    yanchor = "bottom",
                    xanchor = "left",
                    orientation = "h",
                ),
            ),
        )
    return p1, p2 
end 

function _plot_box_plot_mean_errors(dataset_to_compare, results_to_compare; selected_points = nothing)
    traces = GenericTrace{Dict{Symbol, Any}}[]
for (ix, r) in enumerate(results_to_compare)
    file = joinpath(r.exp_folder, "input_data", string("train_", r.train_id, ".json"))
    params = TrainParams(file)
    if dataset_to_compare == "train"
        dataset_id = params.train_data.id
        data_ground_truth = Serialization.deserialize(params.train_data_path)
    elseif dataset_to_compare == "validation"
        dataset_id = params.validation_data.id
        data_ground_truth = Serialization.deserialize(params.validation_data_path)
    elseif dataset_to_compare == "test"
        dataset_id = params.test_data.id
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
                string(dataset_id),
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

    #Don't include 0.0 values 
    ir_mean_error = filter(x->x != 0.0, ir_mean_error)
    ii_mean_error = filter(x->x != 0.0, ii_mean_error)
    @show length(vcat(ir_mean_error, ii_mean_error))
    @error "mean error" Statistics.mean(vcat(ir_mean_error, ii_mean_error))
    x0_data =
        vcat(repeat(["real current"], length(ir_mean_error)), repeat(["imag. current"], length(ii_mean_error)))
    @error length( vcat(ir_mean_error, ii_mean_error))
    if selected_points === nothing
        trace1 = box(;
        y = vcat(ir_mean_error, ii_mean_error),
        x = x0_data,
        boxpoints = "all",
        name = r.name,
        marker_color = r.color,
        legendgroup = "1",
        marker_size = 2,
    )
    else 
        trace1 = box(;
        y = vcat(ir_mean_error, ii_mean_error),
        x = x0_data,
        boxpoints = "all",
        name = r.name,
        marker_color = r.color,
        legendgroup = "1",
        selectedpoints = selected_points[ix],
        selected_marker_color = "black",
        selected_marker_size= 5,
        marker_size = 2,
    )
    end 
    push!(traces, trace1)
end
layout = Layout(;
    font_family = "Times New Roman",
    font_size = 18,
    xaxis = attr(font_size = 12, showline=true, linewidth=1, linecolor="black"),
    yaxis = attr(title = "MAE (p.u. current)",
                 tickvals = [1e-4, 1e-3, 0.01, 0.1], 
                 ticktext=["1e-4", "1e-3", "1e-2", "1e-1"], 
                 font_size = 12, 
                 zeroline = false, 
                 automargin = true,
                 showline=true, 
                 type="log", 
                 linewidth=1, 
                 linecolor="black"
                 ),
    legend = attr(
        x = 0.2,
        y = 1.0,
        font = attr(size=14),
        orientation = "h",
        yanchor = "bottom",
        xanchor = "left",
    ),
    boxmode = "group",
    template = "plotly_white",
    #margin = (b =20, t=0, r = 0, l =50)
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
        dataset_id = params.train_data.id
        data_ground_truth = Serialization.deserialize(params.train_data_path)
    elseif dataset_to_compare == "validation"
        dataset_id = params.validation_data.id
        data_ground_truth = Serialization.deserialize(params.validation_data_path)
    elseif dataset_to_compare == "test"
        dataset_id = params.test_data.id
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
                string(dataset_id),
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
        marker_color = r.color,
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
    name = "Ground Truth",
    legendgroup = "1",
    marker_color = :black,
    marker_size = 3,
    #jitter = 100,
    #pointpos = 0,
)
push!(traces, trace1)

layout = Layout(;
    font_family = "Times New Roman",
    font_size = 16,
    xaxis = attr(font_size = 12, showline=true, automargin =true, linewidth=1, linecolor="black"),
    yaxis = attr(title = "Simulation time (s)", 
    tickvals = [1, 2, 5, 10], 
    #ticktext=["1e-4", "1e-3", "1e-2", "1e-1"], 
    font_size = 12, zeroline = false, type = "log", showline=true, linewidth=1, linecolor="black"),
    legend = attr(
        x = 0.25,
        y = 1.0,
        font = attr(size=14),
        yanchor = "bottom",
        xanchor = "middle",
        orientation = "h",
    ),
    boxmode = "group",
    template = "plotly_white",
    margin = (b =35, t=30, r = 50, l =50)
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
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    for (ix, d) in enumerate(data_ground_truth)
        if d.stable == true 
            p = make_subplots(
                rows = 2,
                cols = 1,
                specs = [Spec(); Spec()][:, 1:1],
               # subplot_titles = ["Real Current (p.u.)" "Imaginary Current (p.u.)" ],
                vertical_spacing = 0.20,

            )
            traces_2 = GenericTrace{Dict{Symbol, Any}}[] 
            for (i,r) in enumerate(results_to_compare)
                file =  joinpath(r.exp_folder, "input_data", string("train_", r.train_id, ".json"))
                params = TrainParams(file)
                dataset_surrogate = Serialization.deserialize(joinpath(params.output_data_path, params.train_id, string("surrogate_", dataset_to_compare, "_", string(params.test_data.id), "_", string(r.chosen_iteration), "_dataset")))
                if  dataset_surrogate[ix].stable ==true 
                    trace_surrogate_ir = PlotlyJS.scatter(;x=dataset_surrogate[ix].tsteps, y = get_device_terminal_data(dataset_surrogate[ix])[:ir],  mode="lines",  name =r.name, line_color = r.color)
                    trace_surrogate_ii = PlotlyJS.scatter(;x=dataset_surrogate[ix].tsteps, y = get_device_terminal_data(dataset_surrogate[ix])[:ii],  mode="lines",  name =r.name, line_color = r.color, showlegend=false)
                    add_trace!(p, trace_surrogate_ir, row =1, col=1)
                    add_trace!(p, trace_surrogate_ii, row =2, col=1)

                end 
            end 
            trace_ground_truth_ir =  PlotlyJS.scatter(;x=d.tsteps, y= get_device_terminal_data(d)[:ir], mode="lines", name = "Ground Truth", line_color = "black")
            trace_ground_truth_ii =  PlotlyJS.scatter(;x=d.tsteps, y= get_device_terminal_data(d)[:ii], mode="lines", name = "Ground Truth", line_color = "black", showlegend=false)
            add_trace!(p, trace_ground_truth_ir, row =1, col=1)
            add_trace!(p, trace_ground_truth_ii, row =2, col=1)


            relayout!(p, showlegend = true)
            p.plot.layout.xaxis = attr(automargin=true,  font_size = 12,  zeroline = false,  linecolor="black" )
            p.plot.layout.xaxis2 = attr(title = "Time (s)", font_size =12, linecolor="black")
            p.plot.layout.yaxis =  attr(automargin = true, title = "Real current (p.u.)", font_size =12, linecolor ="black", zeroline = false) # range = [-1.85, -1.5]
            p.plot.layout.yaxis2 = attr(title = "Imag. current (p.u.)", font_size =12, linecolor="black",  zeroline = false)
            p.plot.layout.template = "plotly_white"
            p.plot.layout.font_family = "Times New Roman"
            p.plot.layout.font_size = 16
            p.plot.layout.legend = attr(
                x = 0.2,
                y = 1.05,
                font_size = 14,
                yanchor = "bottom",
                xanchor = "middle",
                orientation = "h",
            )
            #display(p)
            if ix == 3  #Arbitrary choice of trace to include in paper
                PlotlyJS.savefig(
                    p,
                    joinpath(@__DIR__, "..", "outputs", string("traces", ix, ".pdf")), 
                    width = 400,
                    height = 450,
                    scale = 1,
                )
                PlotlyJS.savefig(
                    p,
                    joinpath(@__DIR__, "..", "outputs", string("traces", ix, ".png")), 
                    width = 400,
                    height = 450,
                    scale = 1,
                )
            end 
        end 
    end  
end 

