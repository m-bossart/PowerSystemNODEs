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
##
#const 
train_split = 0.9999
optimizer = ADAM(0.01)#BFGS() #Optim.KrylovTrustRegion() #ADAM(0.01) #BFGS()
solver = Rodas4() #KenCarp4() #QBDF()# TRBDF2() #Rodas4() #Rodas5() TRBDF2() 
abstol = 1e-6
reltol = 1e-3
tfault =  0.01
tspan = (0.0, 1.0)
steps = 50
tsteps =  10 .^ (range(log10(tfault), log10(tspan[2]),length= steps))
group_size = 50 #5
batching_factor = 1
scale_maxmin = 10
lb_loss = 0.005
nn_width = 2
maxiters = 500
nn_hidden = 1 
nn_activation = gelu  #tanh #
nn_scale = 1.0  #1e-1, 1e-2
n_checkpoint = 10 
is_restart = false
display_plots = false  




label = "run1"
nn_hidden = 2 
nn_width = 3
include("train_nn.jl")

label = "run2"
nn_hidden = 2 
nn_width = 4
include("train_nn.jl")

label = "run3"
nn_hidden = 2 
nn_width = 5
include("train_nn.jl")

label = "run4"
nn_hidden = 2 
nn_width = 6
include("train_nn.jl")

label = "run5"
nn_hidden = 3 
nn_width = 3
include("train_nn.jl")

label = "run6"
nn_hidden = 3 
nn_width = 4
include("train_nn.jl")

label = "run7"
nn_hidden = 3 
nn_width = 5
include("train_nn.jl")

label = "run8"
nn_hidden = 3 
nn_width = 6
include("train_nn.jl")

label = "run9"
nn_hidden = 4 
nn_width = 3
include("train_nn.jl")

label = "run10"
nn_hidden = 4 
nn_width = 4
include("train_nn.jl")

label = "run11"
nn_hidden = 4 
nn_width = 5
include("train_nn.jl")

label = "run12"
nn_hidden = 4 
nn_width = 6
include("train_nn.jl")
