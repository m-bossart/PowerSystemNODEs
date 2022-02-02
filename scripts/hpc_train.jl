using PowerSimulationNODE
using Base.Cartesian

function build_params_list!(params_data, no_change_params, change_params)
    train_id = 1
    starting_dict = no_change_params
    dims = []
    for (k, v) in change_params
        push!(dims, length(v))
    end
    dims_tuple = tuple(dims...)
    iterator = CartesianIndices(dims_tuple)
    for i in iterator
        for (j, (key, value)) in enumerate(change_params)
            starting_dict[key] = value[i[j]]
        end
        starting_dict[:train_id] = string(train_id)
        push!(params_data, NODETrainParams(; starting_dict...))
        starting_dict = no_change_params
        train_id += 1
    end
end

params_data = NODETrainParams[]
no_change_params = Dict{Symbol, Any}()
change_params = Dict{Symbol, Any}()

#INDICATE CONSTANT, NON-DEFAULT PARAMETERS
no_change_params[:maxiters] = 2000
no_change_params[:node_width] = 20
no_change_params[:node_layers] = 2
no_change_params[:graphical_report_mode] = 3

#INDICATE PARAMETES TO ITERATE OVER COMBINATORIALLY 
change_params[:node_activation] = ["relu", "hardtanh", "sigmoid"]
change_params[:training_groups] = [
    [(
        tspan = (0.0, 1.0),
        multiple_shoot_group_size = 10,
        multiple_shoot_continuity_term = 100,
        batching_sample_factor = 1.0,
    )],
    [(
        tspan = (0.0, 1.0),
        multiple_shoot_group_size = 20,
        multiple_shoot_continuity_term = 100,
        batching_sample_factor = 1.0,
    )],
    [(
        tspan = (0.0, 1.0),
        multiple_shoot_group_size = 101,
        multiple_shoot_continuity_term = 100,
        batching_sample_factor = 1.0,
    )],
]
change_params[:node_state_inputs] = [
    [],
    [("gen1", :ir_filter), ("gen1", :ii_filter)],
    [("gen1", :q_oc), ("gen1", :θ_oc), ("gen1", :ω_oc)],
    [
        ("gen1", :γd_ic),
        ("gen1", :γq_ic),
        ("gen1", :ξd_ic),
        ("gen1", :ξq_ic),
        ("gen1", :ϕd_ic),
        ("gen1", :ϕq_ic),
    ],
    [
        ("gen1", :vi_filter),
        ("gen1", :γd_ic),
        ("gen1", :vq_pll),
        ("gen1", :γq_ic),
        ("gen1", :ξd_ic),
        ("gen1", :ϕd_ic),
        ("gen1", :ε_pll),
        ("gen1", :ir_cnv),
        ("gen1", :vr_filter),
        ("gen1", :ω_oc),
        ("gen1", :ξq_ic),
        ("gen1", :vd_pll),
        ("gen1", :q_oc),
        ("gen1", :ϕq_ic),
        ("gen1", :θ_pll),
        ("gen1", :θ_oc),
        ("gen1", :ii_cnv),
    ],
]

build_params_list!(params_data, no_change_params, change_params)
@warn length(params_data)
##
#=
 hpc_params = SavioHPCTrain(;
    username = "jdlara",
    params_data = params_data,
    project_folder = "PowerSystemNODEs",
    scratch_path = "/global/home/users/jdlara",
)
  =#

hpc_params = SummitHPCTrain(;
    username = "mabo4366",
    params_data = params_data,
    project_folder = "PowerSystemNODEs",
    scratch_path = "/scratch/summit/mabo4366",
    n_tasks = length(params_data),
    force_generate_inputs = true,
)

generate_train_files(hpc_params)
run_parallel_train(hpc_params)
