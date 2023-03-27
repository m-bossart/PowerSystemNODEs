using PowerSystems
const PSY = PowerSystems
using PowerFlows
using PowerSimulationsDynamics
#using Sundials
using Random

function _compute_total_load_parameters(load::PSY.StandardLoad)
    # Constant Power Data
    constant_active_power = PSY.get_constant_active_power(load)
    constant_reactive_power = PSY.get_constant_reactive_power(load)
    max_constant_active_power = PSY.get_max_constant_active_power(load)
    max_constant_reactive_power = PSY.get_max_constant_reactive_power(load)
    # Constant Current Data
    current_active_power = PSY.get_current_active_power(load)
    current_reactive_power = PSY.get_current_reactive_power(load)
    max_current_active_power = PSY.get_max_current_active_power(load)
    max_current_reactive_power = PSY.get_max_current_reactive_power(load)
    # Constant Admittance Data
    impedance_active_power = PSY.get_impedance_active_power(load)
    impedance_reactive_power = PSY.get_impedance_reactive_power(load)
    max_impedance_active_power = PSY.get_max_impedance_active_power(load)
    max_impedance_reactive_power = PSY.get_max_impedance_reactive_power(load)
    # Total Load Calculations
    active_power = constant_active_power + current_active_power + impedance_active_power
    reactive_power =
        constant_reactive_power + current_reactive_power + impedance_reactive_power
    max_active_power =
        max_constant_active_power + max_current_active_power + max_impedance_active_power
    max_reactive_power =
        max_constant_reactive_power +
        max_current_reactive_power +
        max_impedance_reactive_power
    return active_power, reactive_power, max_active_power, max_reactive_power
end

function transform_load_to_constant_impedance(load::PSY.StandardLoad)
    # Total Load Calculations
    active_power, reactive_power, max_active_power, max_reactive_power =
        _compute_total_load_parameters(load)
    # Set Impedance Power
    PSY.set_impedance_active_power!(load, active_power)
    PSY.set_impedance_reactive_power!(load, reactive_power)
    PSY.set_max_impedance_active_power!(load, max_active_power)
    PSY.set_max_impedance_reactive_power!(load, max_reactive_power)
    # Set everything else to zero
    PSY.set_constant_active_power!(load, 0.0)
    PSY.set_constant_reactive_power!(load, 0.0)
    PSY.set_max_constant_active_power!(load, 0.0)
    PSY.set_max_constant_reactive_power!(load, 0.0)
    PSY.set_current_active_power!(load, 0.0)
    PSY.set_current_reactive_power!(load, 0.0)
    PSY.set_max_current_active_power!(load, 0.0)
    PSY.set_max_current_reactive_power!(load, 0.0)
    return
end

sys = System("systems/36Bus_CR.json")
for l in get_components(StandardLoad, sys)
    transform_load_to_constant_impedance(l)
end 

power_flow_results_pre = run_powerflow(sys)

line = get_component(DynamicBranch, sys, "Bus 6-Bus 26-i_1")
set_x!(line.branch, get_x(line.branch) * 3.5)
set_b!(line.branch, (from = get_b(line.branch).from * 3, to = get_b(line.branch).to * 3))
new_line = deepcopy(line.branch)
set_name!(new_line, "Bus 6-Bus 26-i_2")
add_component!(sys, new_line)
power_flow_results_post = run_powerflow(sys)

line = get_component(DynamicBranch, sys, "Bus 5-Bus 15-i_1")
set_x!(line.branch, get_x(line.branch) * 2.5)

line = get_component(DynamicBranch, sys, "Bus 38-Bus 39-i_1")
set_x!(line.branch, get_x(line.branch) * 4)

line = get_component(DynamicBranch, sys, "Bus 35-Bus 34-i_1")
set_x!(line.branch, get_x(line.branch) * 4)

load = get_component(StandardLoad, sys, "load181")
set_impedance_active_power!(load, get_impedance_active_power(load) * 1.5)

load = get_component(StandardLoad, sys, "load281")
set_impedance_active_power!(load, get_impedance_active_power(load) * 1.2)

gen = get_component(ThermalStandard, sys, "generator-3-Trip")
set_active_power!(gen, get_active_power(gen) + 0.3)

load = get_component(StandardLoad, sys, "load361")
set_impedance_active_power!(load, get_impedance_active_power(load) * 0.8)

load = get_component(StandardLoad, sys, "load61")
set_impedance_active_power!(load, get_impedance_active_power(load) * 0.75)

power_flow_results_post = run_powerflow(sys)

gfm_bats = get_components(
    x -> isa(get_freq_estimator(get_dynamic_injector(x)), KauraPLL),
    GenericBattery,
    sys,
)

for b in gfm_bats
    @show gfm_available = round(rand()) > 0
    gfl_bat = collect(get_components(x -> get_bus(x) == get_bus(b), GenericBattery, sys))
    gfl_bat =
        filter(x -> !isa(get_freq_estimator(get_dynamic_injector(x)), KauraPLL), gfl_bat)
    if gfm_available
        set_active_power!(b, get_active_power(b) * 2)
    else
        set_available!(b, gfm_available)
        set_active_power!(first(gfl_bat), get_active_power(b) * 2)
    end
end

set_active_power!(get_component(GenericBattery, sys, "GFM_Battery-26"), 1.0)

set_active_power!(get_component(GenericBattery, sys, "Gfl_Battery-22"), 1.0)

set_magnitude!(get_bus(get_component(GenericBattery, sys, "Gfl_Battery-21")), 1.018)

run_powerflow!(sys)

gf_bats = get_components(
    x -> isa(get_freq_estimator(get_dynamic_injector(x)), KauraPLL),
    GenericBattery,
    sys,
)

for b in gf_bats
    !get_available(b) && continue
    bus = get_bus(b)
    set_bustype!(bus, BusTypes.PQ)
    @show get_name(b)
    @show get_reactive_power(b)
    @show get_active_power(b)
    @show cos(atan(get_reactive_power(b) / get_active_power(b)))
    #set_reactive_power!(b, 0.0)
end

gfm_bats = get_components(
    x -> !isa(get_freq_estimator(get_dynamic_injector(x)), KauraPLL),
    GenericBattery,
    sys,
)

for b in gfm_bats
    !get_available(b) && continue
    bus = get_bus(b)
    set_bustype!(bus, BusTypes.PV)
    @show get_name(b)
    @show get_reactive_power(b)
    @show get_active_power(b)
    @show cos(atan(get_reactive_power(b) / get_active_power(b)))
    #set_reactive_power!(b, 0.0)
end

run_powerflow!(sys)

sim = Simulation(ResidualModel, sys, mktempdir(), (0.0, 10.0))
sm = small_signal_analysis(sim)
summary_eigenvalues(sm)

to_json(sys, "systems/36Bus.json"; force = true)
