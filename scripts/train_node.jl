include("../src/PowerSystemNODEs.jl")
include("../system_data/dynamic_components_data.jl")
configure_logging(console_level = Logging.Info)
#configure_logging(;filename = "train_node.log")

sample_train_parameters = "input_data/sample_train_parameters.json"
serialize(NODETrainParams(), "input_data/sample_train_parameters.json")

train_params_file = isempty(ARGS) ? sample_train_parameters : ARGS[1]
train_params = NODETrainParams(train_params_file)
train_params.graphical_report = true
status = train(train_params)    #compare to previous serial version 

#LOCAL TEST OF MULTIPLE RUNS BELOW FOR TESTING summarize_trains.jl
train_params.train_id = "train_instance_2"
train_params.maxiters = 6
status = train(train_params)
