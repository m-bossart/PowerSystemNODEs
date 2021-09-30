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

train_split = 0.9999
optimizer = ADAM(0.01)  #BFGS(), Optim.KrylovTrustRegion(),
solver = Rodas4()       #KenCarp4(), QBDF(), TRBDF2() 
abstol = 1e-6
reltol = 1e-3
tfault =  0.1
tspan = (0.0, 2.0)
steps = 150
tsteps =  10 .^ (range(log10(tfault), log10(tspan[2]),length= steps))  
tsteps = tspan[1]:((tspan[2]-tspan[1])/steps):tspan[2]
group_size = 150 
batching_factor = 1 #NOT IMPLEMENTED YET 
lb_loss = 0.000000005
nn_width = 5
maxiters = 5
nn_hidden = 1 
nn_activation = gelu  #tanh 
nn_scale = 1.0  #1e-1, 1e-2
plot_log = false 
display_plots = false   
loss_function = "mae"  # "mae" 

##
label = "run1"
loss_function = "mae"
nn_hidden = 5 
nn_width = 5
group_size = 150
maxiters = 2000
optimizer = ADAM(0.01)
include("train_nn.jl")

label = "run2"
loss_function = "mae"
nn_hidden = 5 
nn_width = 5
group_size = 75
maxiters = 1000
optimizer = ADAM(0.01)
include("train_nn.jl")

label = "run3"
loss_function = "mae"
nn_hidden = 5 
nn_width = 5
group_size = 25
maxiters = 333
optimizer = ADAM(0.01)
include("train_nn.jl")

##


