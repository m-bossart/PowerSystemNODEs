include("../src/PowerSystemNODEs.jl")
include("../system_data/dynamic_components_data.jl")
configure_logging(console_level = Logging.Info)
#configure_logging(;filename = "train_node.log")

sample_train_parameters = "input_data/sample_train_parameters.json"
serialize(NODETrainParams(), "input_data/sample_train_parameters.json")

train_params_file = isempty(ARGS) ? sample_train_parameters : ARGS[1]
train_params = NODETrainParams(train_params_file)

status = train(train_params)

if train_params.graphical_report
    plots = visualize_training(train_params_1)
    #TODO save plots, and move inside train() function 
end
