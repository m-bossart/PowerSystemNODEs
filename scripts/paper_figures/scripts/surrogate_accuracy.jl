#The main result of the paper. 
using PlotlyJS
using DataFrames
using Plots
using PowerSimulationsDynamics
using PowerSystems
using PowerSimulationNODE
using Serialization
using LaTeXStrings
using JSON3

include(joinpath(@__DIR__, "surrogate_accuracy_plot_utils.jl"))
mkpath(joinpath(@__DIR__, "..", "outputs"))
########### INPUT DATA ###########
dataset_to_compare = "test"
results_to_compare = [
    (
        exp_folder = "data_from_hpc/05_22_23_data_grid",
        train_id = "002",
        chosen_iteration = 0,
        name = "data-driven surrogate-02",
        generate_data = false,
    ),
    (
        exp_folder = "data_from_hpc/05_22_23_data_grid",
        train_id = "003",
        chosen_iteration = 0,
        name = "data-driven surrogate-03",
        generate_data = false,
    ),
    (
        exp_folder = "data_from_hpc/05_22_23_data_grid",
        train_id = "004",
        chosen_iteration = 0,
        name = "data-driven surrogate-4",
        generate_data = false,
    ),
#=      (
            exp_folder = "transfers/exp_04_16_23_physics_grid",
            train_id = "007",
            chosen_iteration = 0,
            name = "physics-based surrogate (trained)",
            generate_data = false,
        ),
        (
            exp_folder = "transfers/exp_04_16_23_physics_grid",
            train_id = "007",
            chosen_iteration = 1,
            name = "physics-based surrogate (untrained)",
            generate_data = false,
        ), =#
]

#Fix all the paths if the results have been copied to a new directory
for r in results_to_compare
    file =  joinpath(pwd(), r.exp_folder, PowerSimulationNODE.INPUT_FOLDER_NAME, string("train_", r.train_id, ".json"))
    rebase_path!(file, r.exp_folder)
end 

# REGENERATE DATASETS 
_regenerate_datasets(dataset_to_compare, results_to_compare)

# PLOT HISTOGRAM OF ERRORS
p_ir, p_ii = _plot_historgram_all_errors(dataset_to_compare, results_to_compare)
#display(p_ir)
PlotlyJS.savefig(
    p_ir, 
    joinpath(@__DIR__, "..", "outputs", "histogram_ir.png"),
    width = 500,
    height = 400,
    scale = 1,
)
#display(p_ii)
PlotlyJS.savefig(
    p_ii, 
    joinpath(@__DIR__, "..", "outputs", "histogram_ii.png"),
    width = 500,
    height = 400,
    scale = 1,
)

# PLOT BOX PLOTS OF MEAN ERRORS 
p_accuracy = _plot_box_plot_mean_errors(dataset_to_compare, results_to_compare)
#display(p_accuracy)
PlotlyJS.savefig(
    p_accuracy,
    joinpath(@__DIR__, "..", "outputs", "box_plot_node_error.png"),
    width = 500,
    height = 400,
    scale = 1,
)


p_timing = _plot_timing_comparison(dataset_to_compare, results_to_compare)
#display(p_timing)
PlotlyJS.savefig(
    p_timing,
    joinpath(@__DIR__, "..", "outputs", "box_plot_times.png"),
    width = 500,
    height = 400,
    scale = 1,
)


# LOOK AT INDIVIDUAL DATA TRACES 
#_display_comparisons_individual_traces(dataset_to_compare, results_to_compare)




##
#Parameter change plots 
#= for r in results_to_compare
    file = joinpath(r.exp_folder, "input_data", string("train_", r.train_id, ".json"))
    show_parameter_change(file)
end =#

