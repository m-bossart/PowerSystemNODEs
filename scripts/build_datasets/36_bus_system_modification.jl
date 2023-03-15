using PowerSystems
using PowerFlows
using PowerSimulationsDynamics
using Sundials
using Random

sys = System("systems/36Bus.json")

power_flow_results_pre = run_powerflow(sys)

line = get_component(DynamicBranch, sys, "Bus 6-Bus 26-i_1")
set_x!(line.branch, get_x(line.branch)*3.5)
set_b!(line.branch, (from = get_b(line.branch).from*3, to = get_b(line.branch).to*3))
new_line = deepcopy(line.branch)
set_name!(new_line, "Bus 6-Bus 26-i_2")
add_component!(sys, new_line)
power_flow_results_post = run_powerflow(sys)

line = get_component(DynamicBranch, sys, "Bus 5-Bus 15-i_1")
set_x!(line.branch, get_x(line.branch)*2.5)

line = get_component(DynamicBranch, sys, "Bus 38-Bus 39-i_1")
set_x!(line.branch, get_x(line.branch)*4)

line = get_component(DynamicBranch, sys, "Bus 35-Bus 34-i_1")
set_x!(line.branch, get_x(line.branch)*4)

load = get_component(PowerLoad, sys, "load181")
set_active_power!(load, get_active_power(load)*1.5)

load = get_component(PowerLoad, sys, "load281")
set_active_power!(load, get_active_power(load)*1.2)

gen = get_component(ThermalStandard, sys, "generator-3-Trip")
set_active_power!(gen, get_active_power(gen) + 0.3)

load = get_component(PowerLoad, sys, "load361")
set_active_power!(load, get_active_power(load)*0.8)

load = get_component(PowerLoad, sys, "load61")
set_active_power!(load, get_active_power(load)*0.75)

power_flow_results_post = run_powerflow(sys)

gfm_bats = get_components(x -> isa(get_freq_estimator(get_dynamic_injector(x)), KauraPLL), GenericBattery, sys)

for b in gfm_bats
    @show gfm_available = round(rand()) > 0
    gfl_bat = collect(get_components(x -> get_bus(x) == get_bus(b), GenericBattery, sys))
    gfl_bat = filter(x -> !isa(get_freq_estimator(get_dynamic_injector(x)), KauraPLL), gfl_bat)
    if gfm_available
        set_active_power!(b, get_active_power(b)*2)
    else
        set_available!(b, gfm_available)
        set_active_power!(first(gfl_bat), get_active_power(b)*2)
    end
end

set_active_power!(get_component(GenericBattery, sys, "GF_Battery-26"), 1.0)

set_active_power!(get_component(GenericBattery, sys, "Gf_Battery-22"), 1.0)

sim = Simulation(ResidualModel, sys, mktempdir(), (0.0, 10.0))
sm = small_signal_analysis(sim)
summary_eigenvalues(sm)

to_json(sys, "systems/36Bus.json"; force=true)
