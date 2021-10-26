
mutable struct NODETrainParams  #add train id. Name the folder based on id 
    solver::String
    solver_tols::Tuple{Float64, Float64}
    sensealg::String
    optimizer::String
    optimizer_η::Float64
    optimizer_adjust::String
    optimizer_adjust_η::Float64
    maxiters::Int64   #should this be for the full training or per call to the optimizer
    lb_loss::Float64  #default to 0 
    batching::Bool
    batch_factor::Float64 #multiply the size of the current training range by batch_factor to get the number of points.
    batch_seed::Int64     #re-create results?  #ONLY NEED ONE SEED 
    groupsize_steps::Int64  #the size of the extending timespan         
    groupsize_faults::Int64 #how many faults should be trained on simultaneously, default to 1. 
    loss_function_weights::Tuple{Float64, Float64}
    loss_function_scale::String     #"none", "range" 
    ode_model::String
    node_input_scale::Float64   #scale on the voltage input only
    node_output_scale::Float64  #scale on the current output only 
    node_inputs::String
    node_feedback_states::Int64
    node_feedback_current::Bool
    node_layers::Int64
    node_width::Int64
    node_activation::String
    node_seed::Int64
    export_mode::Int64
    export_file_path::String
end

StructTypes.StructType(::Type{NODETrainParams}) = StructTypes.Struct()

function NODETrainParams(;
    solver = "Rodas4",
    solver_tols = (1e-6, 1e-9),
    sensealg = "ForwardDiffSensitivity",
    optimizer = "Adam",
    optimizer_η = 0.01,
    optimizer_adjust = "nothing",
    optimizer_adjust_η = 0.01,
    maxiters = 5,
    lb_loss = 0.0,
    batching = false,
    batch_factor = 1.0,
    batch_seed = 1234,
    groupsize_steps = 55,
    groupsize_faults = 1,
    loss_function_weights = (0.5, 0.5),
    loss_function_scale = "range",
    ode_model = "vsm",
    node_input_scale = 1e1,
    node_output_scale = 1.0,
    node_inputs = "voltage",
    node_feedback_states = 0,
    node_feedback_current = true,
    node_layers = 2,
    node_width = 2,
    node_activation = "relu",
    node_seed = 1234,
    export_mode = 3,
    export_file_path = joinpath(pwd(), "data"),
)

    #HERE IS THE LOGIC OF FILLING IN SOME OF THE PARAMETERS THAT MIGHT NOT MAKE SENSE       
    NODETrainParams(
        solver,
        solver_tols,
        sensealg,
        optimizer,
        optimizer_η,
        optimizer_adjust,
        optimizer_adjust_η,
        maxiters,
        lb_loss,
        batching,
        batch_factor,
        batch_seed,
        groupsize_steps,
        groupsize_faults,
        loss_function_weights,
        loss_function_scale,
        ode_model,
        node_input_scale,
        node_output_scale,
        node_inputs,
        node_feedback_states,
        node_feedback_current,
        node_layers,
        node_width,
        node_activation,
        node_seed,
        export_mode,
        export_file_path,
    )
end

struct TrainData
    sys::System
    faults_results_dir::String
end

function read_data(pvs, df)
    df_row = df[df.id .== get_name(pvs), :]     #add a warning in case data doesn't match 
    id = df_row[1, :id]
    tsteps = df_row[1, :tsteps]
    i_true = vcat(df_row[1, :ir_true]', df_row[1, :ii_true]')
    i_ver = vcat(df_row[1, :ir_ver]', df_row[1, :ii_ver]')
    p_ode = df_row[1, :p_ode]
    x₀ = df_row[1, :x₀]
    p_V₀ = df_row[1, :V₀]
    return id, tsteps, i_true, i_ver, p_ode, x₀, p_V₀
end

function train(params::NODETrainParams, data::TrainData)
    #INSTANTIATE 
    sensealg = instantiate_sensealg(params)
    solver = instantiate_solver(params)
    optimizer = instantiate_optimizer(params)
    nn = instantiate_nn(params)
    M = instantiate_M(params)
    !(params.optimizer_adjust == "nothing") &&
        (optimizer_adjust = instantiate_optimizer_adjust(params))

    #READ DATA AND INPUT SYSTEM 
    df = DataFrame(Arrow.Table(data.faults_results_dir))    #change to dictionary 
    pvss = collect(get_components(PeriodicVariableSource, data.sys))

    #TRAIN SEQUENTIALLY 
    local min_θ, res, output_data #res = nothing - remove local 
    output_data = []
    min_θ = initial_params(nn)
    for pvs in pvss
        @show min_θ[1]
        res, output_data =
            train(min_θ, params, sensealg, solver, optimizer, nn, M, df, pvs, output_data)
        min_θ = copy(res.u)
    end

    #TRAIN ADJUSTMENTS (TO DO)

    (params.export_mode == 1) && push!(output_data, (min_θ, res.minimum))   #export level = 1, only save final parameters and loss. 

    capture_output(output_data, params.export_file_path, params.export_mode)

    return res
end

function train(θ, params, sensealg, solver, optimizer, nn, M, df, pvs, output_data) #move 159-168
    #READ FAULT DATA FOR THE CURRENT PVS 
    id, tsteps, i_true, i_ver, p_ode, x₀, p_V₀ = read_data(pvs, df)

    #DEFINE FUNCTIONS OF TIME FOR THE CURRENT PVS
    Vm, Vθ = Source_to_function_of_time(pvs)
    surr = instantiate_surr(params, nn, Vm, Vθ)

    #DEFINE SCALE FACTORS FOR THE CURRENT PVS , REMOVE 
    Vr_scale =
        1 / (
            maximum(Vm.(tsteps) .* cos.(Vθ.(tsteps))) -
            minimum(Vm.(tsteps) .* cos.(Vθ.(tsteps)))
        )
    Vi_scale =
        1 / (
            maximum(Vm.(tsteps) .* sin.(Vθ.(tsteps))) -
            minimum(Vm.(tsteps) .* sin.(Vθ.(tsteps)))
        )

    #FIND INITIAL CONDITIONS FOR THE SURROGATE, BUILD AND SOLVE THE PROBLEM, AND CONFIRM IT MATCHES THE VERIFICATION DATA PROVIDED. ABSTRACT TO FUNCITON (170-191), call inside train 
    p_scale = [Vr_scale, Vi_scale, 0.0]     #turn off the nn 
    p_nn = initial_params(nn)
    n_weights_nn = length(p_nn)
    p_fixed = vcat(p_ode, p_scale, p_V₀, n_weights_nn)
    p = vcat(p_nn, p_fixed)
    order_surr = ODE_ORDER + 2 + params.node_feedback_states + 2   #2 states for the node current, 2 algebraic states   
    x₀_surr = zeros(order_surr)
    x₀_surr[1:ODE_ORDER] = x₀
    x₀_surr[(end - 1):end] = [x₀[I__IR_FILTER], x₀[I__II_FILTER]]
    h = get_init_vsm_v_t_0(p, x₀[I__IR_FILTER], x₀[I__II_FILTER], nn, Vm, Vθ) #make generic! get rid of hard code 
    res_surr = nlsolve(h, x₀_surr)
    @assert converged(res_surr)
    dx = similar(x₀_surr)
    surr(dx, res_surr.zero, p, 0.0)
    @assert all(isapprox.(dx, 0.0; atol = 1e-8))

    tspan = (tsteps[1], tsteps[end])   #tspan comes from data   #Pass flag as parameter in NODETrainParams if you want to verify or not 
    surr_func = ODEFunction(surr, mass_matrix = M)
    surr_prob = ODEProblem(surr_func, x₀_surr, tspan, p)
    sol = solve(
        surr_prob,
        solver,
        abstol = params.solver_tols[1],
        reltol = params.solver_tols[2],
        saveat = tsteps,
    )
    @assert mae(sol[22, :], i_ver[1, :]) < 5e-5

    #PREPARE THE SURROGATE FOR TRAINING  
    p_scale = [Vr_scale, Vi_scale, params.node_output_scale]
    p_fixed = vcat(p_ode, p_scale, p_V₀, n_weights_nn)
    p = vcat(p_nn, p_fixed)
    surr_prob = remake(surr_prob, p = p)

    #Calculate scale factors for loss function based on true data - Create a function for all of the scaling. 
    Ir_scale = maximum(i_true[1, :]) - minimum(i_true[1, :])
    Ii_scale = maximum(i_true[2, :]) - minimum(i_true[2, :])
    #@show surr_init = instantiate_surr_init(params) #don't need? just use surr? 

    u₀ = res_surr.zero

    #INSTANTIATE THE TRAINING FUNCTIONS 
    pred_function = instantiate_pred_function(
        p_fixed,
        solver,
        surr_prob,
        params.solver_tols,
        sensealg,
        u₀,
    )
    loss_function = instantiate_loss_function(
        params.loss_function_weights,
        Ir_scale,
        Ii_scale,
        pred_function,
    )

    try
        #TRAIN ON A SINGLE FAULT USING EXTENDING TIME RANGES.
        datasize = length(tsteps)
        ranges = extending_ranges(datasize, params.groupsize_steps)
        local min_θ, res #remove local 
        min_θ = θ     #Copy? Equal? 
        range_count = 1
        for range in ranges
            @show min_θ[1]
            i_curr = i_true[:, range]
            t_curr = tsteps[range]
            train_loader = Flux.Data.DataLoader(
                (i_curr, t_curr),
                batchsize = Int(floor(length(i_curr[1, :]))),
            )     #TODO - IMPLEMENT BATCHING
            optfun = OptimizationFunction(
                (θ, p, batch, time_batch) -> loss_function(θ, batch, time_batch),
                GalacticOptim.AutoForwardDiff(),
            )
            optprob = OptimizationProblem(optfun, min_θ)

            cb = instantiate_cb!(
                output_data,
                params.lb_loss,
                params.export_mode,
                id,
                range_count,
            )
            range_count += 1

            res = GalacticOptim.solve(
                optprob,
                optimizer,
                ncycle(train_loader, params.maxiters),
                cb = cb,
            )
            min_θ = copy(res.u)
            @show min_θ[1]
            @assert res.minimum == loss_function(res.u, i_curr, t_curr)[1]
            @assert res.minimum == loss_function(min_θ, i_curr, t_curr)[1]
        end
        @show min_θ
        return res, output_data
    catch e
        return 0
    end
end

function capture_output(inputs, export_file_path, export_mode)
    df = DataFrame(inputs)

    (export_mode == 3) && rename!(df, [:id, :range_count, :l, :p, :pred_1, :pred_2])
    (export_mode == 2) && rename!(df, [:id, :range_count, :l])
    (export_mode == 1) && rename!(df, [:p, :l])

    open(joinpath(export_file_path, "outputdata"), "w") do io
        Arrow.write(io, df)
    end
end
