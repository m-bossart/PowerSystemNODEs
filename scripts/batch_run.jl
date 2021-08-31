########################PARAMETERS##############################################
 #base label for training figures
#TODO Improve the plotting of power system results (PowerGraphics?)
using Pkg
Pkg.activate(".")
using Revise
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
using GalacticOptim

#const 
train_split = 0.9999
optimizer = ADAM(0.1)
solver = Rodas4() #KenCarp4() #QBDF()# TRBDF2() #Rodas4() #Rodas5() TRBDF2() 
abstol = 1e-6
reltol = 1e-3
tfault =  0.01
tspan = (0.0, 1.0)
steps = 100
tsteps =  10 .^ (range(log10(tfault), log10(tspan[2]),length= steps))
group_size = 10 
lb_loss = 0.05 
nn_width = 2
maxiters = 500
nn_hidden = 1 
#CHANGE NN DEPTH MANUALLY (default = 1 hidden)! - build a function that can build 1-5 depth nns! 
nn_activation = relu  #tanh #gelu
nn_scale = 1.0  #1e-1, 1e-2

label = "hidden=1,width=3,group=10"
nn_hidden = 1 
nn_width = 3    #do NOT use nn_width<3, see: https://github.com/m-bossart/PowerSystemUDEs/issues/11
include("train_nn.jl")

label = "hidden=1,width=4,group=10"
nn_hidden = 1 
nn_width = 4
include("train_nn.jl")

label = "hidden=1,width=5,group=10"
nn_hidden = 1 
nn_width = 5
include("train_nn.jl")


label = "hidden=2,width=3,group=10"
nn_hidden = 2 
nn_width = 3
include("train_nn.jl")

label = "hidden=2,width=4,group=10"
nn_hidden = 2 
nn_width = 4
include("train_nn.jl")

label = "hidden=2,width=5,group=10"
nn_hidden = 2 
nn_width = 5
include("train_nn.jl")