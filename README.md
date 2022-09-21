# PowerSystemUDEs

Repository of scripts for training SteadyStateNODE surrogates.

## Workflow

At the top of each script, have a "training directory" user input --> might make sense to have global variable that you can set, then you only need to set it once?? 
The training should create the folder and all the inputs/outputs should be stored there in a self-contained folder (except the system? - make the names descriptive)

```
scripts
    analyze 
        summarize_trains.jl
        visualize_trains.jl 
    archive (old scripts)
    build_datasets (scripts for building specific datasets and examining properties independent of training)
    build_systems
        build_9bus.jl
        build_14bus.jl
    hpc_train
        hpc_train.jl (the script that you run to do the training)
        generate_data.jl (called during hpc_train.jl)
        train_node.jl (called during hpc_train.jl)
        build_subsystems.jl (called during hpc_train.jl)
    local_train
        localone script with all of the commands (can easily jump back and forth)
system_data
    dynamic_components_data.jl
systems
    9bus.json
    14bus.json 
```


To get detailed information on the options available use: `? NODETrainParams`

File structure in the training directory:
```
input_data
   data.json
   system_validation_descriptors.json
   system.json   
output_data
    <train_id_1>
        output.json
        output.arrow
    <train_id_2>
        output.json
        output.arrow
    ...
train_parameters
    <train_id_1>.json
    <train_id_2>.json
    ...
```

To run training on hpc:
1) Modify `scripts/hpc_train/hpc_train.jl` with parameters.
2) Run `julia --project scripts/hpc_train/hpc_train.jl` from the PowerSystemNODEs directory. 

Note: to run on hpc, set environment variable  `ENV["GKSwstype"] = "100"` in startup file: `~/.julia/config/startup.jl`