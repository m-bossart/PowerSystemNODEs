using PowerSystems
using PowerSimulationsDynamics
using PowerSimulationNODE
include(joinpath(pwd(), "system_data", "dynamic_components_data.jl"))
raw_file_path =
    joinpath(PowerSimulationNODE.INPUT_SYSTEM_FOLDER_NAME, "IEEE 14 bus_modified_33.raw")
base_system_path = joinpath("systems", "IEEE_14bus_modified.json")

sys = node_load_system(raw_file_path)
surrogate_bus_number = 16
surrogate_bus =
    collect(get_components(Bus, sys, x -> get_number(x) == surrogate_bus_number))[1]
gens = collect(
    get_components(
        ThermalStandard,
        sys,
        x -> get_number(get_bus(x)) == surrogate_bus_number,
    ),
)
!(length(gens) == 1) && @error "number of devices at surrogate bus not equal to one"
gen = gens[1]
total_rating = get_rating(gen) #doesn't impact dynamics
total_base_power = get_base_power(gen)
total_active_power = get_active_power(gen)
remove_component!(sys, gen)
n_inverters = 1
for i in 1:n_inverters
    g = ThermalStandard(
        name = string("gen", string(i)),
        available = true,
        status = true,
        bus = surrogate_bus,
        active_power = total_active_power, #Only divide base power by n_devices
        reactive_power = 0.0,
        rating = total_rating / n_inverters,
        active_power_limits = (min = 0.0, max = 3.0),
        reactive_power_limits = (-3.0, 3.0),
        ramp_limits = nothing,
        operation_cost = ThreePartCost(nothing),
        base_power = total_base_power / n_inverters,
    )
    add_component!(sys, g)
    if (i == 1)
        inv_typ = inv_case78(get_name(g))
        add_component!(sys, inv_typ, g)
    end
    if (i == 2)
        inv_typ = inv_gfoll(get_name(g))
        add_component!(sys, inv_typ, g)
    end
    if (i == 3)
        inv_typ = inv_darco_droop(get_name(g))
        add_component!(sys, inv_typ, g)
    end
end

#Add dynamic models to non-surrogate generators 
gens = collect(
    get_components(
        ThermalStandard,
        sys,
        x -> get_number(get_bus(x)) != surrogate_bus_number,
    ),
)
for g in gens
    inv_typ = inv_case78(get_name(g))
    add_component!(sys, inv_typ, g)
end

to_json(sys, base_system_path, force = true)
