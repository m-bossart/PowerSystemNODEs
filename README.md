# PowerSystemNODEs

Repository for paper: "Acceleration of Power System Dynamic Simulations using a Deep Equilibrium Layer and Neural ODE Surrogate" (https://arxiv.org/abs/2405.06827)  

This repository requires Julia 1.7.

The `paper_experiments`folder contains the scripts for generating the results in the paper on the CU Boulder hpc. 
The `paper_figures` includes the analysis and plotting routines for generating the figures. 

To get detailed information on the training options available use: `? NODETrainParams`
To run an experiment run: `julia --project paper_experiments/<experiment_script.jl>` from the PowerSystemNODEs directory.

Note: experiments for hyper-parameter tuning are designed to be run in parallel in an hpc environment. 
Note: to run on hpc, set environment variable  `ENV["GKSwstype"] = "100"` in startup file: `~/.julia/config/startup.jl` (`/projects/mabo4366/.julia/config`)