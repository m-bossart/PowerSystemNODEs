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

#const 
train_split = 0.9999
optimizer = ADAM(0.01)#BFGS() #Optim.KrylovTrustRegion() #ADAM(0.01) #BFGS()
solver = Rodas4() #KenCarp4() #QBDF()# TRBDF2() #Rodas4() #Rodas5() TRBDF2() 
abstol = 1e-6
reltol = 1e-3
tfault =  0.1#0.01#0.1
tspan = (0.0, 2.0)
steps = 150
tsteps =  10 .^ (range(log10(tfault), log10(tspan[2]),length= steps))  
tsteps = tspan[1]:((tspan[2]-tspan[1])/steps):tspan[2]
#tsteps = 0.001:((tspan[2]-tspan[1])/steps):tspan[2]
#tsteps = tfault:((tspan[2]-tfault)/steps):tspan[2]
#tsteps = vcat(0.0,tsteps) #add 0 for better plots, check ss 
group_size = 200 #5
batching_factor = 1
scale_maxmin = 10
lb_loss = 0.00005
nn_width = 2
maxiters = 2
nn_hidden = 1 
nn_activation = gelu  #tanh #
nn_scale = 1.0  #1e-1, 1e-2

display_plots = true  
loss_function = "mse"  # "mse" 

#What runs do we need for paper?
#Group size: 5, 20, 50, 200 (keep maxiters constant)
#nn_width = 3, 5, 7 
#nn_depth = 3, 5, 7 

label = "test1"
loss_function = "mse"
nn_hidden = 5 
nn_width = 5
group_size = 250
maxiters = 5
optimizer = ADAM(0.01)
include("train_nn.jl")



