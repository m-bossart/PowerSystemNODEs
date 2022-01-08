include("../src/PowerSystemNODEs.jl")
include("../src/constants.jl")
include("../system_data/dynamic_components_data.jl")
configure_logging(console_level = Logging.Info)
#configure_logging(;filename = "train_node.log")

sample_train_parameters = "input_data/sample_parameters.json"
p = NODETrainParams()
p.verify_psid_node_off = false
serialize(p, "input_data/sample_parameters.json")

train_params_file = isempty(ARGS) ? sample_train_parameters : ARGS[1]
train_params = NODETrainParams(train_params_file)
status = train(train_params)
