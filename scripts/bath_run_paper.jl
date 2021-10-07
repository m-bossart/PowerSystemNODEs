########################PARAMETERS##############################################
 #base label for training figures
#TODO Improve the plotting of power system results (PowerGraphics?)
using Pkg
Pkg.activate(".")
using Revise
using BSON: @save, @load
using Distributions
using OrdinaryDiffEq
using DifferentialEquations
using DiffEqSensitivity
using PowerSystems
using Logging
using PowerSimulationsDynamics
const PSID = PowerSimulationsDynamics
const PSY = PowerSystems
using Plots
using FFTW
using IterTools
using Statistics
using NLsolve
using DiffEqFlux: group_ranges
using DiffEqFlux
using Flux
using Flux.Losses: mae, mse
using ForwardDiff
using Statistics
using Optim
using GalacticOptim
using DelimitedFiles

include("../models/constants.jl")
include("../models/DynamicComponents.jl")
include("../models/SurrogateModels.jl")
include("../models/utils.jl")
include("../models/parameter_utils.jl")
include("../models/init_functions.jl")
configure_logging(console_level = Logging.Error)



train_split = 0.9999
optimizer = ADAM(0.01)  #BFGS(), Optim.KrylovTrustRegion(),
solver = Rodas4()       #KenCarp4(), QBDF(), TRBDF2() 
abstol = 1e-6
reltol = 1e-3
tfault =  0.1
tspan = (0.0, 2.0)
steps = 150
group_size = 150 
tsteps =  10 .^ (range(log10(tfault), log10(tspan[2]),length= steps))  
tsteps = tspan[1]:((tspan[2]-tspan[1])/steps):tspan[2]
lb_loss = 0
nn_width = 5
nn_hidden = 5 
maxiters = 1000  #Change 
nn_activation = gelu  
nn_scale = 1.0 
plot_log = false 
display_plots = false   
loss_function = "mae" 

#Voltage is the only external input to the NN, Increase the # of feedback states 
label = "nn_v_2"
surr = vsm_nn_v_2
surr_init = get_init_vsm_nn_v_2
n_extra = 2 #number of nn outputs 
M = MassMatrix(21, 2, 0)
nn = build_nn(4, 2, nn_width, nn_hidden, nn_activation)
include("train_nn.jl")

label = "nn_v_3"
surr = vsm_nn_v_3
surr_init = get_init_vsm_nn_v_3
n_extra = 3 
M = MassMatrix(21, 2, 1)
nn = build_nn(5, 3, nn_width, nn_hidden, nn_activation)
include("train_nn.jl")

label = "nn_v_4"
surr = vsm_nn_v_4
surr_init = get_init_vsm_nn_v_4
n_extra = 4 
M = MassMatrix(21, 2, 2)
nn = build_nn(6, 4, nn_width, nn_hidden, nn_activation)
include("train_nn.jl")

label = "nn_v_5"
surr = vsm_nn_v_5
surr_init = get_init_vsm_nn_v_5
n_extra = 5
M = MassMatrix(21, 2, 3)
nn = build_nn(7, 5, nn_width, nn_hidden, nn_activation) 
include("train_nn.jl")


#= 
#Voltage and Iode current are external input to the NN, Increase the # of feedback states 
label = "nn_vi_2"
surr = vsm_nn_vi_2
M = MassMatrix(21, 2, 0)
nn = build_nn(6, 2, nn_width, nn_hidden, nn_activation)
include("train_nn.jl")

label = "nn_v_3"
surr = vsm_nn_vi_3
M = MassMatrix(21, 2, 1)
nn = build_nn(7, 3, nn_width, nn_hidden, nn_activation) 
include("train_nn.jl")

label = "nn_v_4"
surr = vsm_nn_vi_2
M = MassMatrix(21, 2, 2)
nn = build_nn(8, 4, nn_width, nn_hidden, nn_activation) 
include("train_nn.jl")

label = "nn_v_5"
surr = vsm_nn_vi_2
M = MassMatrix(21, 2, 3)
nn = build_nn(9, 5, nn_width, nn_hidden, nn_activation)  
include("train_nn.jl")
 =#
 
