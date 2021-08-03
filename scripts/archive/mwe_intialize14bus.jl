using Pkg
Pkg.activate(".")
using Revise
using OrdinaryDiffEq
using PowerSystems
using PowerSimulationsDynamics
include("../models/DynamicComponents.jl")

sys = System("cases/IEEE 14 bus_modified_33.raw")
sys2 = System("cases/IEEE 14 bus_modified_33_RemoveFixedAdmittance.raw")

Ybus1 = Ybus(sys)
Ybus2 = Ybus(sys2)
test = Ybus1.data - Ybus2.data
for g in get_components(ThermalStandard,sys)
    inv_case = inv_case78(get_name(g))
    add_component!(sys, inv_case, g)
end

sim = Simulation!(MassMatrixModel, sys, pwd(), (0.0,1.0))
display(solve_powerflow(sys)["bus_results"])
display(get_initial_conditions(sim)["Vm"])

for b in get_components(Bus,sys2)
    @info get_number(b), get_magnitude(b)
end

for g in get_components(ThermalStandard,sys2)
    inv_case = inv_case78(get_name(g))
    add_component!(sys2, inv_case, g)
end

sim = Simulation!(MassMatrixModel, sys2, pwd(), (0.0,1.0))
display(solve_powerflow(sys2)["bus_results"])
display(get_initial_conditions(sim)["Vm"])

for b in get_components(Bus,sys2)
    @info get_number(b), get_magnitude(b)
end
