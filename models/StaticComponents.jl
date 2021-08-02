StaticGen(name::String, bus::Bus) = ThermalStandard(
    name = name,
    available = true,
    status = true,
    bus = bus,
    active_power = 0.0,
    reactive_power = 0.0,
    rating = 0.0,
    active_power_limits=(min=0.0, max=0.0),
    reactive_power_limits=nothing,
    ramp_limits=nothing,
    operation_cost=ThreePartCost(nothing),
    base_power=100.0,
)
