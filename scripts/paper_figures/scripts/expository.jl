#https://stackoverflow.com/questions/69016568/unable-to-export-plotly-images-to-png-with-kaleido
using Pkg
exp_folder = "data_from_hpc/exp_data_grid_expository"
#Pkg.activate(exp_folder)
#Pkg.instantiate()
using PowerSimulationNODE
using Serialization
using LaTeXStrings
using PlotlyJS

########### INPUT DATA ########### 

train_id = "010"
pre_train_iteration = 1
post_train_iteration = 950
dataset_iteration = 1
y_scale_states = (-1.5, 1.5)
y_scale_ir = (-1.9, -1.5)
y_scale_ii = (0, 0.2)
################################# 
#_regenerate_datasets(exp_folder)

df_predictions = PowerSimulationNODE.read_arrow_file_to_dataframe(
    joinpath(exp_folder, "output_data", train_id, "predictions"),
)
surrogate_prediction_before = df_predictions[pre_train_iteration, "surrogate_solution"]
surrogate_prediction_after = df_predictions[post_train_iteration, "surrogate_solution"]
for (i, surrogate_prediction) in
    enumerate([surrogate_prediction_before, surrogate_prediction_after])
    train_dataset =
        Serialization.deserialize(joinpath(exp_folder, "input_data", "train_data_1"))[dataset_iteration]
    @error fieldnames(typeof(train_dataset))
    @error train_dataset.device_terminal_data
    tsteps = train_dataset.tsteps
    ground_truth_real_current = train_dataset.device_terminal_data["Bus 6 -> Bus 26"][:ir]
    ground_truth_imag_current = train_dataset.device_terminal_data["Bus 6 -> Bus 26"][:ii]

    #lay = Plots.@layout [c{0.5w} [a; b]]

    r0_pred = surrogate_prediction.r0_pred
    t_series = surrogate_prediction.t_series
    i_series = surrogate_prediction.i_series
    i_series = reshape(i_series, (2, length(t_series)))       #Loses shape when serialized/deserialized to arrow  
    r_series = surrogate_prediction.r_series
    dim_r = Int(length(r_series) / length(t_series))
    r_series = reshape(r_series, (dim_r, length(t_series)))    #Loses shape when serialized/deserialized to arrow

    p = make_subplots(
        rows = 3,
        cols = 1,
        specs = [Spec(); Spec(); Spec()][:, 1:1], #weird indexing makes matrix not vector    ,
        #subplot_titles = ["Neural ODE States" "Real Current (p.u.)" "Imaginary Current (p.u.)" missing],
        vertical_spacing = 0.08,
        horizontal_spacing = 0.15,
    )
    add_trace!(
        p,
        PlotlyJS.scatter(;
            x = t_series,
            y = r_series[1, :],
            mode = "lines",
            :line => Dict(:color => "gray", :width => 2),
            showlegend = true,
            name = "neural ODE states",
        ),
        row = 1,
        col = 1,
    )
    add_trace!(
        p,
        PlotlyJS.scatter(;
            x = t_series,
            y = ground_truth_real_current,
            showlegend = true,
            :line => Dict(:color => "black", :width => 3, :dash => "solid"),
            mode = "lines",
            name = "ground truth output",
            xaxis = "x1",
            yaxis = "y1",
        ),
        row = 2,
        col = 1,
    )
    add_trace!(
        p,
        PlotlyJS.scatter(;
            x = zeros(dim_r),
            y = r0_pred,
            mode = "markers",
            :line => Dict(:color => "#2ca02c", :width => 3),
            showlegend = true,
            name = "predicted initial conditions",
        ),
        row = 1,
        col = 1,
    )
    add_trace!(
        p,
        PlotlyJS.scatter(;
            x = t_series,
            y = i_series[1, :],
            mode = "lines",
            showlegend = true,
            name = "surrogate predicted output",
            :line => Dict(:color => "#2ca02c", :dash => "dashdot", :width => 3),
        ),
        row = 2,
        col = 1,
    )

    add_trace!(
        p,
        PlotlyJS.scatter(;
            x = t_series,
            y = ground_truth_imag_current,
            mode = "lines",
            showlegend = false,
            :line => Dict(:color => "black", :dash => "solid", :width => 3),
        ),
        row = 3,
        col = 1,
    )

    add_trace!(
        p,
        PlotlyJS.scatter(;
            x = t_series,
            y = i_series[2, :],
            mode = "lines",
            showlegend = false,
            :line => Dict(:color => "#2ca02c", :width => 3, :dash => "dashdot"),
            xaxis = "x2",
            yaxis = "y2",
        ),
        row = 3,
        col = 1,
    )

    for i in 1:size(r_series, 1)
        add_trace!(
            p,
            PlotlyJS.scatter(;
                x = t_series,
                y = r_series[i, :],
                mode = "lines",
                :line => Dict(:color => "gray", :width => 2),
                showlegend = false,
                xaxis = "x3",
                yaxis = "y3",
            ),
            row = 1,
            col = 1,
        )
    end


    relayout!(p, showlegend = true)
    p.plot.layout.xaxis = attr( font_size = 12, autpomargin=true,  zeroline = false,  linecolor="black" )
    p.plot.layout.xaxis2 = attr( zeroline = false ,  font_size = 12, linecolor="black")
    p.plot.layout.xaxis3 = attr(title = "Time (s)",  font_size =12,  zeroline = false , linecolor="black")
    p.plot.layout.yaxis =  attr(title = "Neural ODE States", font_size =12, linecolor ="black", zeroline = false)
    p.plot.layout.yaxis2 = attr(title = "Real current (p.u.)", font_size =12, linecolor="black")
    p.plot.layout.yaxis3 = attr(title = "Imag. current (p.u.)", font_size =12, linecolor="black", automargin=true ,zeroline = false)
    p.plot.layout.template = "plotly_white"
    p.plot.layout.font_family = "Times New Roman"
    p.plot.layout.font_size = 18
    p.plot.layout.legend = attr(
        x = -0.15,
        y = 1.05,
        font_size = 17,
        yanchor = "bottom",
        xanchor = "left",
        orientation = "h",
    )
    #p.plot.layout.xaxis_domain=[0, 10.0],
    #yaxis_domain=[0, 0.45],
    #xaxis4=attr(domain=[0.55, 1.0], anchor="y4"),
    #xaxis2_domain=[0.55, 1],
    #yaxis3_domain=[0.55, 1],
    #yaxis4=attr(domain=[0.55, 1], anchor="x4")
    if i == 1
        p.plot.layout.showlegend = true 
        PlotlyJS.savefig(
            p,
            joinpath(@__DIR__, "..", "outputs", "expository_before.pdf"),
            scale = 1.0,
            height = 800,
            width = 500,
        )
        PlotlyJS.savefig(
            p,
            joinpath(@__DIR__, "..", "outputs", "expository_before.png"),
            scale = 1.0,
            height = 800,
            width = 500,
        )
    else 
        p.plot.layout.showlegend = true
        #p.plot.layout.margin = (t=0)
        PlotlyJS.savefig(
            p,
            joinpath(@__DIR__, "..", "outputs", "expository_after.pdf"),
            scale = 1.0,
            height = 800,
            width = 500,
        )
        PlotlyJS.savefig(
            p,
            joinpath(@__DIR__, "..", "outputs", "expository_after.png"),
            scale = 1.0,
            height = 800,
            width = 500,
        )
    end

    #=     p1 = PlotlyJS.plot([p1_trace1, p1_trace2], Layout(;title="Real current", legend =true) )
        p2 = PlotlyJS.plot([p2_trace1, p2_trace2])
        p3 = PlotlyJS.plot(p3_traces)
        display([p3; p1 p2]) =#
end
