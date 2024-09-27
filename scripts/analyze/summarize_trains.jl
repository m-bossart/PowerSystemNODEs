#This can be used to find the best training and iteration across an experiment. 
using PowerSimulationNODE
using Plots

#train_folder = joinpath("transfers", "exp_08_25_23_physics_grid") 
train_folder = joinpath("transfers", "exp_08_23_23_data_random")
train_folder = joinpath("data_from_hpc", "exp_data_sample_0.1")
train_folder = joinpath("data_from_hpc", "exp_data_sample_0.05")
train_folder = joinpath("data_from_hpc", "exp_data_sample_0.2")


plotlyjs()

a = generate_summary(joinpath(pwd(), train_folder, "output_data"))
print_train_parameter_overview(
    joinpath(pwd(), train_folder, PowerSimulationNODE.INPUT_FOLDER_NAME),
)
p = visualize_summary(a)


#Check all validation losses, not just the final loss. 
#Build vector of tuples (train_id, iteration, validation_loss) so we can select the lowest overall -> compare to screenshot to see if we get relatively better 
using Statistics
list_validation_results = []
for dir in readdir(joinpath(pwd(), train_folder, "output_data"), join=true)
    df_validation_loss = PowerSimulationNODE.read_arrow_file_to_dataframe(joinpath(dir, "validation_loss"))
    #print(df_validation_loss)
    for (i, row) in enumerate( eachrow( df_validation_loss ) ) 
        if 0.0 ∉ row[:mae_ir] && 0.0 ∉ row[:mae_ii]
            ir_mean_over_faults = Statistics.mean(row[:mae_ir]) 
            ii_mean_over_faults = Statistics.mean(row[:mae_ii]) 
            push!(list_validation_results, (dir, i, (ir_mean_over_faults + ii_mean_over_faults)/2))
        else 
        end 
    end 
end 

minimum([x[3] for x in list_validation_results ])
min, index = findmin([x[3] for x in list_validation_results ])
list_validation_results[index]

display(list_validation_results[index])
##
p = TrainParams(joinpath(pwd(), "transfers", "exp_08_23_23_data_random", "input_data", "train_185.json" ))
p.optimizer[1].log_η
p.optimizer[1].loss_function.α
n_states = p.model_params.dynamic_hidden_states
p.model_params.initializer_n_layer
3 + p.model_params.initializer_width_layers_relative_input      #width hidden layers
p.model_params.dynamic_n_layer
n_states + p.model_params.dynamic_width_layers_relative_input   #width hidden layers
#TODO- find total number of parameters. 