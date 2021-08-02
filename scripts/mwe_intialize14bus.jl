using Pkg
Pkg.activate(".")
using Revise
using Distributions
using OrdinaryDiffEq
using PowerSystems
using PowerSimulationsDynamics
const PSID = PowerSimulationsDynamics
const PSY = PowerSystems
using Plots
using FFTW

include("../models/DynamicComponents.jl")
include("../models/InverterModels.jl")
include("../models/StaticComponents.jl")
include("../models/utils.jl")
include("../models/init_functions.jl")


sys = System("cases/IEEE 14 bus_modified_33.raw")
tspan = (0.0,1.0)
for g in get_components(ThermalStandard,sys)
    inv_case = inv_case78(get_name(g))
    add_component!(sys, inv_case, g)
end
solve_powerflow(sys)["bus_results"]
sim = Simulation!(MassMatrixModel, sys, pwd(), tspan)
print_device_states(sim)
collect(get_components(DynamicInverter,sys))[1]

x₀_dict = get_initial_conditions(sim)
