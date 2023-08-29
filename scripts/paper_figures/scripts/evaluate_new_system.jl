#= This script should take in a new 36bus system which has been modified in some way,
and show how well the physics based surrogate will perform without any training,
only by predicting a base-power weighted average of parameters. The goal is to find 
a system where simply defining an aggregate surrogate doesn't perform very well and 
we can show the benefit of a data-driven approach. 
 =#
using Revise
using PowerFlows
using PowerSystems
using PowerSimulationsDynamics
using PowerSimulationNODE
using PowerSimulationsDynamicsSurrogates
const PSIDS = PowerSimulationsDynamicsSurrogates
using Logging
using Serialization
using Plots
include(joinpath(pwd(), "scripts", "build_datasets", "utils.jl"))
include(joinpath(pwd(), "scripts", "hpc_train", "utils.jl"))
include(joinpath(pwd(), "scripts", "paper_figures", "scripts", "surrogate_accuracy_plot_utils.jl"))

train_folder = "evaluate_new_system"
system_name = "36bus_fix"
project_folder = "PowerSystemNODEs"
scratch_path = joinpath(pwd(), "..")

_copy_full_system_to_train_directory(
    scratch_path,
    project_folder,
    train_folder,
    system_name,
)

p = TrainParams(
    train_id = "BASE",
    surrogate_buses = vcat(21:29, 31:39),
    test_data = (
        id = "1",
        operating_points = repeat(
            [
                RandomOperatingPointXiao(
                    generator_voltage_range = (0.94, 1.06),
                    generator_power_range = (0.0, 1.0),
                    load_multiplier_range = (0.5, 1.5),
                ),
            ],
            3,#50
        ),
         perturbations = repeat(
            [[PSIDS.RandomLoadChange(time = 1.0, load_multiplier_range = (0.0, 2.0))]],
            1,
        ), 
        #perturbations = repeat([[BranchTrip(1.0, Line, "Bus 6-Bus 26-i_2")]],1),
        #perturbations = repeat([[PSIDS.RandomBranchTrip(time = 1.0)]],1),
        params = PSIDS.GenerateDataParams(
            solver = "Rodas5",
            solver_tols = (reltol = 1e-3, abstol = 1e-6),
            tspan = (0.0, 10.0),
            tstops = 0.0:0.1:10.0,
            tsave = 0.0:0.1:10.0,
            formulation = "MassMatrix",
            all_branches_dynamic = false,
            all_lines_dynamic = false,
            seed = 33,
        ),
    ),
    model_params = MultiDeviceParams(name = "source_1"),
    rng_seed = 11,
    output_mode_skip = 1,
    train_time_limit_seconds = 1e9,
    base_path = joinpath(scratch_path, project_folder, train_folder),
    system_path = joinpath(
        scratch_path,
        project_folder,
        train_folder,
        PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME,
        string(system_name, ".json"),
    ),
)

######################################################################################
################################# BUILD AND GENERATE #################################
######################################################################################
build_subsystems(p)
mkpath(joinpath(p.base_path, PowerSimulationNODE.INPUT_FOLDER_NAME))
generate_test_data(p)
test_dataset = Serialization.deserialize(p.test_data_path)
sys_validation = PowerSimulationNODE.node_load_system(p.surrogate_system_path)
sys_validation_aux = deepcopy(sys_validation)
PowerSimulationNODE.add_surrogate_psid!(sys_validation, p.model_params, test_dataset)
to_json(sys_validation, p.modified_surrogate_system_path, force = true)
sys = System(p.system_path)
p_start =  determine_p_start(sys, p.surrogate_buses)
validation_sys = PowerSimulationNODE.node_load_system(p.modified_surrogate_system_path)
validation_sys_aux = PowerSimulationNODE.node_load_system(p.surrogate_system_path)
test_dataset = Serialization.deserialize(p.test_data_path)
data_collection_location = Serialization.deserialize(p.data_collection_location_path)[2]
data_ground_truth = Serialization.deserialize(p.test_data_path)
 data_surrogate = generate_surrogate_dataset(
    validation_sys,
    validation_sys_aux,
    p_start,
    test_dataset,
    p.test_data,
    data_collection_location,
    p.model_params,
)  
  
using PlotlyJS
traces = GenericTrace{Dict{Symbol, Any}}[]
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
)
push!(traces, trace1)
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
    boxmode = "group",
    template = "plotly_white",
    margin = (b =20, t=0, r = 50, l =50)
)

p = PlotlyJS.plot(traces, layout);
PlotlyJS.savefig(
    p, 
    joinpath(pwd(), "scripts", "paper_figures", "outputs", "original_system.png"),
    width = 500,
    height = 400,
    scale = 1,
)


for (ix, d) in enumerate(data_ground_truth)
    if d.stable == true 
        terminal_data_ground_truth = d.device_terminal_data["Bus 6 -> Bus 26"]
        terminal_data_surrogate = data_surrogate[ix].device_terminal_data["Bus 6 -> Bus 26"]
        traces_ir = GenericTrace{Dict{Symbol, Any}}[] 
        trace_surrogate = PlotlyJS.scatter(;x=data_surrogate[ix].tsteps, y=terminal_data_surrogate[:ir], mode="lines", name = string("ir-", ix))
        push!(traces_ir, trace_surrogate)
        trace_ground_truth =  PlotlyJS.scatter(;x=d.tsteps, y=terminal_data_ground_truth[:ir], mode="lines", name = "ground truth", line_color = "black")
        push!(traces_ir, trace_ground_truth)
        display(PlotlyJS.plot(traces_ir))
        traces_ii = GenericTrace{Dict{Symbol, Any}}[] 
        trace_surrogate = PlotlyJS.scatter(;x=data_surrogate[ix].tsteps, y=terminal_data_surrogate[:ii], mode="lines", name = string("ii-", ix))
        push!(traces_ii, trace_surrogate)
        trace_ground_truth =  PlotlyJS.scatter(;x=d.tsteps, y=terminal_data_ground_truth[:ii], mode="lines", name = "ground truth", line_color = "black")
        push!(traces_ii, trace_ground_truth)
        display(PlotlyJS.plot(traces_ii))
    end 
end  

