include("../system_data/dynamic_components_data.jl")
include("../src/PowerSystemNODEs.jl")
configure_logging(console_level = Logging.Info)

force_generate = isempty(ARGS) ? "true" : ARGS[1]
force_generate_inputs = convert(Bool, force_generate == "true")

train_data_path = joinpath(INPUT_FOLDER_NAME, "data.json")
train_system_path = joinpath(INPUT_FOLDER_NAME, "system.json")
full_system_path = joinpath(INPUT_SYSTEM_FOLDER_NAME, "full_system.json")

if (!(isfile(full_system_path)) || force_generate_inputs)
    @warn "Rebuilding full system "
    include("build_full_system.jl")
    force_generate_inputs = true
end

if (!(isfile(train_system_path)) || !(isfile(train_data_path)) || force_generate_inputs)
    @warn "Rebuilding train system and train input data"
    sys_full = System(full_system_path)
    pvs_data = fault_data_generator("scripts/config.yml")
    sys_pvs = build_pvs(pvs_data)
    label_area!(sys_full, [16], "surrogate")
    @assert check_single_connecting_line_condition(sys_full)
    sys_surr = remove_area(sys_full, "1")
    sys_train = build_train_system(sys_surr, sys_pvs, "surrogate")
    to_json(sys_train, joinpath(INPUT_FOLDER_NAME, "system.json"), force = true)
    d = generate_train_data(sys_train, NODETrainDataParams())   #BUG - only works for single fault
    serialize(d, joinpath(INPUT_FOLDER_NAME, "data.json"))
end

######### POST TRAIN GENERATE PREDICTION DATA ########
#= sys_rest = remove_area(sys_full, "surrogate")
sys_reduced = build_reduced_system(sys_rest, NODE, "1")
prediction_data = fault_data_generator(sys_reduced) 
final_loss = loss(ground_truth_data, prediction_data)=#
