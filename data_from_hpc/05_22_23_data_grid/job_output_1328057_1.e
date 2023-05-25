┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 208 of the NLsolve function has a value of 21.861724079007853.
│                   Generator = source_1, state = r1.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:95
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:95
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 208 of the NLsolve function has a value of 21.861724079007853.
│                   Generator = source_1, state = r1.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 208 of the NLsolve function has a value of 21.861724079007853.
│                   Generator = source_1, state = r1.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 208 of the NLsolve function has a value of 21.861724079007853.
│                   Generator = source_1, state = r1.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 208 of the NLsolve function has a value of 21.861724079007853.
│                   Generator = source_1, state = r1.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 208 of the NLsolve function has a value of 21.861724079007853.
│                   Generator = source_1, state = r1.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 208 of the NLsolve function has a value of 21.861724079007853.
│                   Generator = source_1, state = r1.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 208 of the NLsolve function has a value of 21.861724079007853.
│                   Generator = source_1, state = r1.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 208 of the NLsolve function has a value of 21.861724079007853.
│                   Generator = source_1, state = r1.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 208 of the NLsolve function has a value of 21.861724079007853.
│                   Generator = source_1, state = r1.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 208 of the NLsolve function has a value of 21.861724079007853.
│                   Generator = source_1, state = r1.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 208 of the NLsolve function has a value of 21.861724079007853.
│                   Generator = source_1, state = r1.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 208 of the NLsolve function has a value of 21.861724079007853.
│                   Generator = source_1, state = r1.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 209 of the NLsolve function has a value of 21.080439167791504.
│                   Generator = source_1, state = r2.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 208 of the NLsolve function has a value of 21.861724079007853.
│                   Generator = source_1, state = r1.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 208 of the NLsolve function has a value of 21.861724079007853.
│                   Generator = source_1, state = r1.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 208 of the NLsolve function has a value of 21.861724079007853.
│                   Generator = source_1, state = r1.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 208 of the NLsolve function has a value of 21.861724079007853.
│                   Generator = source_1, state = r1.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 208 of the NLsolve function has a value of 21.861724079007853.
│                   Generator = source_1, state = r1.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 208 of the NLsolve function has a value of 21.861724079007853.
│                   Generator = source_1, state = r1.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
┌ Warning: Initialization in SteadyStateNODE failed
└ @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/SurrogateComponents/SteadyStateNODE.jl:306
┌ Error: PowerSimulationsDynamics.MassMatrixModel failed to build
│   exception =
│    The initial residual in 208 of the NLsolve function has a value of 21.861724079007853.
│                   Generator = source_1, state = r1.
│                   Error is too large to continue
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] _check_residual(residual::Vector{Float64}, inputs::PowerSimulationsDynamics.SimulationInputs, tolerance::Float64)
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:111
│      [3] refine_initial_condition!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, model::PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.SimCache{typeof(PowerSimulationsDynamics.system_mass_matrix!)}}, jacobian::PowerSimulationsDynamics.JacobianFunctionWrapper{PowerSimulationsDynamics.var"#67#69"{ForwardDiff.JacobianConfig{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12, Tuple{Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}, Vector{ForwardDiff.Dual{ForwardDiff.Tag{PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Float64}, Float64, 12}}}}, PowerSimulationsDynamics.var"#66#68"{PowerSimulationsDynamics.SystemModel{PowerSimulationsDynamics.MassMatrixModel, PowerSimulationsDynamics.JacobianCache{typeof(PowerSimulationsDynamics.system_mass_matrix!), ForwardDiff.Dual{ForwardDiff.Tag{typeof(PowerSimulationsDynamics.system_mass_matrix!), Float64}, Float64, 12}}}}, Int64}, Matrix{Float64}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/nlsolve_wrapper.jl:136
│      [4] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:405 [inlined]
│      [5] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [6] macro expansion
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:404 [inlined]
│      [7] macro expansion
│        @ /projects/mabo4366/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]
│      [8] _build!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:374
│      [9] (::PowerSimulationsDynamics.var"#108#109"{Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:429
│     [10] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [11] with_logger
│        @ ./logging.jl:623 [inlined]
│     [12] #build!#107
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:428 [inlined]
│     [13] PowerSimulationsDynamics.Simulation(::Type{PowerSimulationsDynamics.MassMatrixModel}, system::System, simulation_folder::String, tspan::Tuple{Float64, Float64}, perturbations::Vector{PowerSimulationsDynamics.Perturbation}; kwargs::Base.Pairs{Symbol, Bool, Tuple{Symbol, Symbol}, NamedTuple{(:all_branches_dynamic, :all_lines_dynamic), Tuple{Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:188
│     [14] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:227
│     [15] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [16] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [17] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [18] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [19] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [20] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [21] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [22] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [23] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [24] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:419
┌ Error: Execution failed
│   exception =
│    The Simulation status is BUILD_FAILED. Can not continue, correct your inputs and build the simulation again.
│    Stacktrace:
│      [1] error(s::String)
│        @ Base ./error.jl:33
│      [2] simulation_pre_step!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:447
│      [3] _execute!(sim::PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, solver::OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}; kwargs::Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}})
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:473
│      [4] (::PowerSimulationsDynamics.var"#116#117"{Base.Pairs{Symbol, Any, NTuple{7, Symbol}, NamedTuple{(:abstol, :reltol, :tstops, :save_everystep, :saveat, :reset_simulation, :enable_progress_bar), Tuple{Float64, Float64, Vector{Float64}, Bool, Vector{Float64}, Bool, Bool}}}, PowerSimulationsDynamics.Simulation{PowerSimulationsDynamics.MassMatrixModel}, OrdinaryDiffEq.Rodas5{0, true, Nothing, typeof(OrdinaryDiffEq.DEFAULT_PRECS), Val{:forward}, true, nothing}})()
│        @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:530
│      [5] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│      [6] with_logger
│        @ ./logging.jl:623 [inlined]
│      [7] #execute!#115
│        @ /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:528 [inlined]
│      [8] fill_surrogate_data!(data::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, sys_train::System, psid_perturbations::Vector{PowerSimulationsDynamics.Perturbation}, data_collection::PowerSimulationsDynamicsSurrogates.GenerateDataParams, data_aux::PowerSimulationsDynamicsSurrogates.SteadyStateNODEData, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:262
│      [9] generate_surrogate_data(sys_main::System, sys_aux::System, perturbations::Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, operating_points::Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, data_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEDataParams, data_collection_params::PowerSimulationsDynamicsSurrogates.GenerateDataParams; dataset_aux::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, surrogate_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationsDynamicsSurrogates /projects/mabo4366/.julia/packages/PowerSimulationsDynamicsSurrogates/ArusJ/src/generate_data/Datasets.jl:128
│     [10] generate_surrogate_dataset(sys_main::System, sys_aux::System, θ::Vector{Float32}, groundtruth_dataset::Vector{PowerSimulationsDynamicsSurrogates.SteadyStateNODEData}, data_params::NamedTuple{(:id, :operating_points, :perturbations, :params), Tuple{String, Vector{PowerSimulationsDynamicsSurrogates.SurrogateOperatingPoint}, Vector{Vector{Union{PowerSimulationsDynamics.Perturbation, PowerSimulationsDynamicsSurrogates.SurrogatePerturbation}}}, PowerSimulationsDynamicsSurrogates.GenerateDataParams}}, data_collection_location::Vector{Tuple{String, Symbol}}, model_params::PowerSimulationsDynamicsSurrogates.SteadyStateNODEParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:206
│     [11] train(params::TrainParams)
│        @ PowerSimulationNODE /projects/mabo4366/.julia/packages/PowerSimulationNODE/pTrnX/src/train/train.jl:947
│     [12] (::var"#1#2")()
│        @ Main /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:20
│     [13] with_logstate(f::Function, logstate::Any)
│        @ Base.CoreLogging ./logging.jl:511
│     [14] with_logger(f::Function, logger::MultiLogger)
│        @ Base.CoreLogging ./logging.jl:623
│     [15] top-level scope
│        @ /scratch/alpine/mabo4366/PowerSystemNODEs/scripts/hpc_train/train_node.jl:18
│     [16] include(mod::Module, _path::String)
│        @ Base ./Base.jl:418
│     [17] exec_options(opts::Base.JLOptions)
│        @ Base ./client.jl:292
│     [18] _start()
│        @ Base ./client.jl:495
└ @ PowerSimulationsDynamics /projects/mabo4366/.julia/packages/PowerSimulationsDynamics/1kjnZ/src/base/simulation.jl:532
