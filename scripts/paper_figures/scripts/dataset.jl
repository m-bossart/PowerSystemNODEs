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
exp_folder = "data_from_hpc/05_22_23_data_grid"
train_id = 1

train_data_path = joinpath(exp_folder, "input_data", string("train_data_", train_id))
validation_data_path =
    joinpath(exp_folder, "input_data", string("validation_data_", train_id))
test_data_path = joinpath(exp_folder, "input_data", string("test_data_", train_id))

train_dataset = Serialization.deserialize(train_data_path)
validation_dataset = Serialization.deserialize(validation_data_path)
test_dataset = Serialization.deserialize(test_data_path)


####### PLOT TRAIN, VALIDATION, AND TEST TOGETHER  ################
plot_train = true
plot_validation = false
plot_test = false
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
    horizontal_spacing = 0.11,
)

if plot_train
    for t in train_dataset
        add_trace!(
            p,
            PlotlyJS.scatter(;
                x = t.tsteps,
                y = vec(t.surrogate_real_voltage),
                line_color = "black",
            ),
            row = 1,
            col = 1,
        )
        add_trace!(
            p,
            PlotlyJS.scatter(;
                x = t.tsteps,
                y = vec(t.surrogate_imag_voltage),
                line_color = "black",
            ),
            row = 1,
            col = 2,
        )
        add_trace!(
            p,
            PlotlyJS.scatter(; x = t.tsteps, y = vec(t.real_current), line_color = "black"),
            row = 2,
            col = 1,
        )
        add_trace!(
            p,
            PlotlyJS.scatter(; x = t.tsteps, y = vec(t.imag_current), line_color = "black"),
            row = 2,
            col = 2,
        )
    end
end
relayout!(p, showlegend = false)
p.plot.layout.xaxis = attr(title = "Time (s)", showgrid = true, zeroline = false)
p.plot.layout.xaxis2 = attr(title = "Time (s)", showgrid = true, zeroline = false)
p.plot.layout.xaxis3 = attr(title = "Time (s)", showgrid = true, zeroline = false)
p.plot.layout.xaxis4 = attr(title = "Time (s)", showgrid = true, zeroline = false)
p.plot.layout.yaxis = attr(title = "Real voltage (p.u.)", showgrid = true, zeroline = false)
p.plot.layout.yaxis2 =
    attr(title = "Imag. voltage (p.u.)", showgrid = true, zeroline = false)
p.plot.layout.yaxis3 =
    attr(title = "Real current (p.u.)", showgrid = true, zeroline = false)
p.plot.layout.yaxis4 =
    attr(title = "Imag current (p.u.)", showgrid = true, zeroline = false)
p.plot.layout.font_size = 18
#display(p)
PlotlyJS.savefig(
    p,
    joinpath(@__DIR__, "..", "outputs", "dataset.png"),
    scale = 1.0,
    height = Int64(2.25 * 300),
    width = Int64(3.5 * 300),
)

