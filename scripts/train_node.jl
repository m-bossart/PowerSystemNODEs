using Revise
using PowerSimulationNODE
using Logging
include("../system_data/dynamic_components_data.jl")

configure_logging(console_level = Logging.Info)
#configure_logging(;filename = "train_node.log")

sample_train_parameters = "input_data/sample_parameters.json"
p = NODETrainParams()
p.verify_psid_node_off = false
PowerSimulationNODE.serialize(p, "input_data/sample_parameters.json")

train_params_file = isempty(ARGS) ? sample_train_parameters : ARGS[1]

#train_params = NODETrainParams(train_params_file)
#status = train(train_params)




#FOR DEBUGGING -- DELETE BELOW
train_params = NODETrainParams(
    train_id = "1",
    maxiters = 10,
    training_groups = [(tspan = (0.0, 1.0), multiple_shoot_group_size = 10, multiple_shoot_continuity_term = 100.0, batching_sample_factor = 1.0)],
    node_width = 40,
    node_layers = 5,
    node_activation = "relu",
    graphical_report_mode = 3,
)
@time train(train_params)

