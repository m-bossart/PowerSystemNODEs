using PowerSimulationNODE
using Logging

train_params_file = ARGS[1]
train_params = TrainParams(train_params_file)
build_subsystems(train_params)
