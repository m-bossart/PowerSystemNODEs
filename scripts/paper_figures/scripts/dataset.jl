using Revise
using PowerSimulationNODE
using PowerSimulationsDynamicsSurrogates
const PSIDS = PowerSimulationsDynamicsSurrogates
using Logging
using Serialization
using LaTeXStrings
using PlotlyJS
include(joinpath(pwd(), "scripts", "build_datasets", "utils.jl"))
include(joinpath(pwd(), "scripts", "hpc_train", "utils.jl"))

########### INPUT DATA ###########      
exp_folder = "transfers/exp_08_07_23_data_random" #"evaluate_new_system", "data_from_hpc/"
train_id = 1

train_data_path = joinpath(exp_folder, "input_data", string("train_data_", train_id))
validation_data_path =
    joinpath(exp_folder, "input_data", string("validation_data_", train_id))
test_data_path = joinpath(exp_folder, "input_data", string("test_data_", train_id))

####### CHOOSE DATASET TO PLOT  ################
dataset = Serialization.deserialize(train_data_path)
#dataset = Serialization.deserialize(validation_data_path)
#dataset = Serialization.deserialize(test_data_path)

traces_vr = GenericTrace{Dict{Symbol, Any}}[]
traces_vi = GenericTrace{Dict{Symbol, Any}}[]
traces_ir = GenericTrace{Dict{Symbol, Any}}[]
traces_ii = GenericTrace{Dict{Symbol, Any}}[]
p = make_subplots(
    rows = 2,
    cols = 2,
    specs = [Spec() Spec(); Spec() Spec()],
    #subplot_titles = ["Neural ODE States" "Real Current (p.u.)" "Imaginary Current (p.u.)" missing],
    vertical_spacing = 0.1,
    horizontal_spacing = 0.15,
)
ix_max_up = 20 #16 index of maximum load step up 
ix_max_down = 20 #index of maximum load step down

for (ix, t) in enumerate(dataset)
    if t.stable
        v_color = "gray"
        i_color = "black"
        if ix == ix_max_up
            v_color = "red"
            i_color = "red"
        elseif ix == ix_max_down
            v_color = "red"
            i_color = "red"
        end 
        add_trace!(
            p,
            PlotlyJS.scatter(;
                x = t.tsteps,
                y = t.device_terminal_data["Bus 6 -> Bus 26"][:vr],
                line_color = v_color,
            ),
            row = 1,
            col = 1,
        )
        add_trace!(
            p,
            PlotlyJS.scatter(;
                x = t.tsteps,
                y = t.device_terminal_data["Bus 6 -> Bus 26"][:vi],
                line_color = v_color,
            ),
            row = 1,
            col = 2,
        )
        add_trace!(
            p,
            PlotlyJS.scatter(; x = t.tsteps, y = t.device_terminal_data["Bus 6 -> Bus 26"][:ir], line_color = i_color),
            row = 2,
            col = 1,
        )
        add_trace!(
            p,
            PlotlyJS.scatter(; x = t.tsteps, y =  t.device_terminal_data["Bus 6 -> Bus 26"][:ii], line_color = i_color),
            row = 2,
            col = 2,
        )
    end 
end

relayout!(p, showlegend = false)
p.plot.layout.xaxis = attr( showgrid = true, zeroline = false, linecolor="black")
p.plot.layout.xaxis2 = attr(  showgrid = true, zeroline = false, linecolor="black")
p.plot.layout.xaxis3 = attr(title = "Time (s)", showgrid = true, zeroline = false, linecolor="black")
p.plot.layout.xaxis4 = attr(title = "Time (s)", showgrid = true, zeroline = false, linecolor="black")
p.plot.layout.yaxis = attr(title = "Real voltage (p.u.)", showgrid = true, zeroline = false, linecolor="black")
p.plot.layout.yaxis2 =
    attr(title = "Imag. voltage (p.u.)", showgrid = true, zeroline = false, linecolor="black")
p.plot.layout.yaxis3 =
    attr(title = "Real current (p.u.)", showgrid = true, zeroline = false, linecolor="black")
p.plot.layout.yaxis4 =
    attr(title = "Imag current (p.u.)", showgrid = true, zeroline = false, linecolor="black")
p.plot.layout.font_size = 12
p.plot.layout.font_family = "Times New Roman"
p.plot.layout.template = "plotly_white"
p.plot.layout.margin = (
    t=20,  #top margin
    b=35,
    r=20,
    l=50,
) 
#display(p)
mkpath(joinpath(@__DIR__, "..", "outputs"))

PlotlyJS.savefig(
    p,
    joinpath(@__DIR__, "..", "outputs", "dataset.pdf"),
    scale = 1.0,
    width = 500,
    height = 375,
    #height = Int64(2.25 * 300),
    #width = Int64(3.5 * 300),
)

