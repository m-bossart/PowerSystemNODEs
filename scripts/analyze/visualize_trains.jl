using PowerSystems
using PowerSimulationNODE
using Plots
using JSON3
using Logging
using Serialization
#include("../system_data/dynamic_components_data.jl")
train_folder = joinpath("transfers", "exp_03_25_23_data_grid")

configure_logging(console_level = Logging.Info)
visualize_level = isempty(ARGS) ? 3 : parse(Int64, ARGS[1])

train_files = filter(
    x -> occursin("train_", x) && occursin(".json", x),
    readdir(
        joinpath(pwd(), train_folder, PowerSimulationNODE.INPUT_FOLDER_NAME),
        join = true,
    ),
)

output_folders = readdir(
    joinpath(pwd(), train_folder, PowerSimulationNODE.OUTPUT_FOLDER_NAME),
    join = false,
)

train_files_with_output = filter(
    x ->
        occursin("train_", x) &&
            occursin(".json", x) &&
            TrainParams(x).train_id in output_folders,
    train_files,
)
##
gr()
run_validation = true
a = time()
for file in [train_files_with_output[52]]
    rebase_path!(file, train_folder)
    params = TrainParams(file)
    path_to_output = joinpath(params.output_data_path, params.train_id)
    output_dict =
        JSON3.read(read(joinpath(path_to_output, "high_level_outputs")), Dict{String, Any})
    n_recorded_iterations = length(output_dict["recorded_iterations"])
    visualize_training(file, vcat(1:50, (n_recorded_iterations - 50):n_recorded_iterations))
    #animate_training(file, skip = 100) 
    if run_validation
        validation_sys = node_load_system(params.modified_surrogate_system_path)
        validation_sys_aux = node_load_system(params.surrogate_system_path)
        validation_dataset = Serialization.deserialize(params.validation_data_path)
        data_collection_location =
            Serialization.deserialize(params.data_collection_location_path)[2]

        #Get the parameters (θ) resulting from the training.    
        df_predictions = PowerSimulationNODE.read_arrow_file_to_dataframe(
            joinpath(path_to_output, "predictions"),
        )
        chosen_iteration_index =
            indexin(output_dict["chosen_iteration"], output_dict["recorded_iterations"])[1]
        θ = df_predictions[chosen_iteration_index, "parameters"][1]
        first_sol = df_predictions[1, "surrogate_solution"]
        #        @warn first_sol.deq_iterations  
        #surrogate_dataset_validation = Serialization.deserialize(string(p.output_data_path, "surrogate_validation_dataset"))  #If already generated the data
        surrogate_dataset_validation = generate_surrogate_dataset(
            validation_sys,
            validation_sys_aux,
            θ,
            validation_dataset,
            params.validation_data,
            data_collection_location,
            params.model_params,
        )

        plots_validation_performance =
            visualize_loss(surrogate_dataset_validation, validation_dataset)
        Serialization.serialize(
            joinpath(
                params.output_data_path,
                params.train_id,
                "surrogate_validation_dataset",
            ),
            surrogate_dataset_validation,
        )
        for (i, p) in enumerate(plots_validation_performance)
            Plots.png(p, joinpath(path_to_output, string("validation_dataset_", i)))
        end
    end
    #Cleanup to be able to delete the arrow files: https://github.com/apache/arrow-julia/issues/61
    #GC.gc() 
end

##
for file in train_files_with_output
    rebase_path!(file, train_folder)
    params = TrainParams(file)
    path_to_output = joinpath(params.output_data_path, params.train_id)
    output_dict =
        JSON3.read(read(joinpath(path_to_output, "high_level_outputs")), Dict{String, Any})

    df_predictions = PowerSimulationNODE.read_arrow_file_to_dataframe(
        joinpath(path_to_output, "predictions"),
    )
    first_sol = df_predictions[1, "surrogate_solution"]
    println(split(file, "\\")[end], "     ", first_sol.deq_iterations)
    println("")
end
