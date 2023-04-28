#test if the returned current is correct when a static line is tripped
#https://github.com/NREL-SIIP/PowerSimulationsDynamics.jl/issues/269
using PowerSystems
using PowerSimulationsDynamics
using Sundials
using Plots
sys = System("systems/36bus.json")

pert = BranchTrip(1.0, Line, "Bus 6-Bus 26-i_2")
sim = Simulation(ResidualModel, sys, mktempdir(), (0.0, 10.0), pert)
sm = small_signal_analysis(sim)
summary_eigenvalues(sm)
show_components(sys, Line)

execute!(sim, IDA())

res = read_results(sim)
i1 = get_real_current_branch_flow(res, "Bus 6-Bus 26-i_1")
i2 = get_real_current_branch_flow(res, "Bus 6-Bus 26-i_2")
plot(i1)
plot!(i2)

##
using PowerSystems
exit()
##
