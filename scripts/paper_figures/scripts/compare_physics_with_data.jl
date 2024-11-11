#The main result of the paper. 

using PlotlyJS
using DataFrames
using Plots
using PowerSimulationsDynamics
using PowerSimulationsDynamicsSurrogates
using PowerSystems
using PowerSimulationNODE
using Serialization
using LaTeXStrings
using JSON3

const PSIDS = PowerSimulationsDynamicsSurrogates

include(joinpath(@__DIR__, "surrogate_accuracy_plot_utils.jl"))
mkpath(joinpath(@__DIR__, "..", "outputs"))
########### INPUT DATA ###########
dataset_to_compare = "test"
new_test_dataset = nothing  #defaults to the test_dataset from TrainParams
results_to_compare = [
    (
        exp_folder = "data_from_hpc/exp_08_23_23_data_random",
        train_id = "185",
        chosen_iteration = 8000,
        name = "Data-driven Surrogate",
        color = "#2ca02c",  
        generate_data = false,
    ),
    (
        exp_folder = "data_from_hpc/exp_09_08_24_Load+GFL+GFM",
        train_id = "001",
        chosen_iteration = 0,
        name = "Load + GFL + GFM",
        color = "#1f77b4",
        generate_data = false,
    ),
    (
        exp_folder = "data_from_hpc/exp_09_08_24_GFM",
        train_id = "001",
        chosen_iteration = 1,
        name = "GFM",
        color = "#ff7f0e",
        generate_data = false,
    ), 
]


#Fix all the paths if the results have been copied to a new directory
for r in results_to_compare
    file =  joinpath(pwd(), r.exp_folder, PowerSimulationNODE.INPUT_FOLDER_NAME, string("train_", r.train_id, ".json"))
    rebase_path!(file, r.exp_folder)
end 

#Finds the parameter file for each of the items in results_to_compare.
#Writes in a new value for test_data and also a new path for test_data_path based on the provided id.
#Generates the new ground-truth dataset if generate_data = true.
original_test_data = (
    id = "1",
    operating_points = repeat(
        [
            RandomOperatingPointXiao(
                generator_voltage_range = (0.94, 1.06),
                generator_power_range = (0.0, 1.0),
                load_multiplier_range = (0.5, 1.5),
            ),
        ],
        100,
    ),
     perturbations = repeat(
        [[PSIDS.RandomLoadChange(time = 1.0, load_multiplier_range = (0.0, 2.0))]],
        1,
    ), 
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
)

new_test_data = (
    id = "2",
    operating_points = repeat(
        [
            RandomOperatingPointXiao(
                generator_voltage_range = (0.94, 1.06),
                generator_power_range = (0.0, 1.0),
                load_multiplier_range = (0.5, 1.5),
            ),
        ],
        100,
    ),
     perturbations = repeat(
        [[PSIDS.LineTrip(time = 1.0, line_name = "Bus 6-Bus 26-i_2")]],
        1,
    ), 
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
)
for r in results_to_compare
    file = joinpath(r.exp_folder, "input_data", string("train_", r.train_id, ".json"))
    params = TrainParams(file)   
end 
function change_test_dataset(results_to_compare, new_test_data)
    for r in results_to_compare
        file = joinpath(r.exp_folder, "input_data", string("train_", r.train_id, ".json"))
        params = TrainParams(file)
        params.test_data = new_test_data
        params.test_data_path = string(rstrip(x->x ∈ ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'], params.test_data_path), new_test_data.id)
        PowerSimulationNODE.serialize(params, file)
        if r.generate_data == true 
            PowerSimulationNODE.generate_test_data(params) 
        end 
    end 
end 
change_test_dataset(results_to_compare, original_test_data) #Need this to change the ID! 


# REGENERATE DATASETS 
_regenerate_datasets(dataset_to_compare, results_to_compare)

# PLOT BOX PLOTS OF MEAN ERRORS 
include(joinpath(@__DIR__, "surrogate_accuracy_plot_utils.jl"))
p_accuracy = _plot_box_plot_mean_errors(dataset_to_compare, results_to_compare; selected_points=[[2, 91], [2, 92], [2, 91]]);
PlotlyJS.savefig(
    p_accuracy,
    joinpath(@__DIR__, "..", "outputs", "box_plot_node_error.pdf"),
    width = 400,
    height = 300,
    scale = 1,
);
PlotlyJS.savefig(
    p_accuracy,
    joinpath(@__DIR__, "..", "outputs", "box_plot_node_error.png"),
    width = 400,
    height = 300,
    scale = 1,
);


p_timing = _plot_timing_comparison(dataset_to_compare, results_to_compare)
PlotlyJS.savefig(
    p_timing,
    joinpath(@__DIR__, "..", "outputs", "box_plot_times.pdf"),
    width = 400,
    height = 300,
    scale = 1,
)
p_timing = _plot_timing_comparison(dataset_to_compare, results_to_compare)
PlotlyJS.savefig(
    p_timing,
    joinpath(@__DIR__, "..", "outputs", "box_plot_times.png"),
    width = 400,
    height = 300,
    scale = 1,
)

# LOOK AT INDIVIDUAL DATA TRACES 
_display_comparisons_individual_traces(dataset_to_compare, results_to_compare)


