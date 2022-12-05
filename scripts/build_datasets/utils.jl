
function add_surrogate_A!(sys, connecting_bus_number)
    loads = collect(
        get_components(
            PowerLoad,
            sys,
            x -> get_number(get_bus(x)) == connecting_bus_number,
        ),
    )
    if length(loads) !== 1
        @error "Must specify a bus with exactly one PowerLoad"
        return
    end
    @assert length(loads) == 1
    l = loads[1]
    connecting_bus =
        collect(get_components(Bus, sys, x -> get_number(x) == connecting_bus_number))[1]
    if length(loads) !== 1
        return
    end
    bus_numbers = get_number.(collect(get_components(Bus, sys)))
    new_bus_number = maximum(bus_numbers) + 1
    b = Bus(
        number = new_bus_number,
        name = "Bus_$new_bus_number",
        bustype = BusTypes.PQ,
        angle = 0.0,
        magnitude = 1.0,
        voltage_limits = (0.0, 2.0),
        base_voltage = get_base_voltage(connecting_bus),
    )
    add_component!(sys, b)
    a = Arc(connecting_bus, b)
    add_component!(sys, a)
    li = Line(
        name = "interconnect-line",
        available = true,
        active_power_flow = 0.0,
        reactive_power_flow = 0.0,
        arc = a,
        r = 0.01,
        x = 0.1,
        b = (from = 0.0, to = 0.0),
        rate = 0.0,
        angle_limits = (min = -1.571, max = 1.571),
    )
    add_component!(sys, li)
    set_bus!(l, b)
    return new_bus_number
end

function add_surrogate_B!(sys, connecting_bus_number)
    loads = collect(
        get_components(
            PowerLoad,
            sys,
            x -> get_number(get_bus(x)) == connecting_bus_number,
        ),
    )
    if length(loads) !== 1
        @error "Must specify a bus with exactly one PowerLoad"
        return
    end
    @assert length(loads) == 1
    l = loads[1]
    connecting_bus =
        collect(get_components(Bus, sys, x -> get_number(x) == connecting_bus_number))[1]
    if length(loads) !== 1
        return
    end
    @warn l
    bus_numbers = get_number.(collect(get_components(Bus, sys)))
    new_bus_number = maximum(bus_numbers) + 1
    b = Bus(
        number = new_bus_number,
        name = "Bus_$new_bus_number",
        bustype = BusTypes.PQ,
        angle = 0.0,
        magnitude = 1.0,
        voltage_limits = (0.0, 2.0),
        base_voltage = get_base_voltage(connecting_bus),
    )
    add_component!(sys, b)
    a = Arc(connecting_bus, b)
    add_component!(sys, a)
    li = Line(
        name = "interconnect-line",
        available = true,
        active_power_flow = 0.0,
        reactive_power_flow = 0.0,
        arc = a,
        r = 0.01,
        x = 0.1,
        b = (from = 0.0, to = 0.0),
        rate = 0.0,
        angle_limits = (min = -1.571, max = 1.571),
    )
    add_component!(sys, li)

    @error "Creating battery and removing load "
    remove_component!(sys, l)
    batt = add_battery(
        sys,
        "battery",
        get_name(b),
        l.base_power,
        -l.active_power,
        -l.reactive_power,
    )
    add_component!(sys, batt)
    inv = inv_case78(batt)
    add_component!(sys, inv, batt)
    return new_bus_number
end

function add_surrogate_C!(sys, connecting_bus_number)
    loads = collect(
        get_components(
            PowerLoad,
            sys,
            x -> get_number(get_bus(x)) == connecting_bus_number,
        ),
    )
    if length(loads) !== 1
        @error "Must specify a bus with exactly one PowerLoad"
        return
    end
    @assert length(loads) == 1
    l = loads[1]
    connecting_bus =
        collect(get_components(Bus, sys, x -> get_number(x) == connecting_bus_number))[1]
    if length(loads) !== 1
        return
    end
    @warn l
    bus_numbers = get_number.(collect(get_components(Bus, sys)))
    new_bus_number = maximum(bus_numbers) + 1
    b = Bus(
        number = new_bus_number,
        name = "Bus_$new_bus_number",
        bustype = BusTypes.PQ,
        angle = 0.0,
        magnitude = 1.0,
        voltage_limits = (0.0, 2.0),
        base_voltage = get_base_voltage(connecting_bus),
    )
    add_component!(sys, b)
    a = Arc(connecting_bus, b)
    add_component!(sys, a)
    li = Line(
        name = "interconnect-line",
        available = true,
        active_power_flow = 0.0,
        reactive_power_flow = 0.0,
        arc = a,
        r = 0.01,
        x = 0.1,
        b = (from = 0.0, to = 0.0),
        rate = 0.0,
        angle_limits = (min = -1.571, max = 1.571),
    )
    add_component!(sys, li)

    @error "Creating battery and removing load "
    remove_component!(sys, l)
    batt = add_battery(
        sys,
        "battery",
        get_name(b),
        l.base_power,
        -l.active_power,
        -l.reactive_power,
    )
    add_component!(sys, batt)
    inv = inv_gfoll(batt)
    add_component!(sys, inv, batt)
    return new_bus_number
end

function set_penetration_level!(sys, Gf, GF)
    syncGen = collect(get_components(Generator, sys))
    ibrGen = collect(get_components(GenericBattery, sys))
    allGen = vcat(syncGen, ibrGen)

    ibrBus = unique([gen.bus.number for gen in ibrGen])
    syncBus = unique([gen.bus.number for gen in syncGen])
    genBus = unique(vcat(ibrBus, syncBus))

    busCap = Dict(zip(genBus, zeros(length(genBus))))

    for bus in genBus
        busCap[bus] = sum([
            get_base_power(gen) for gen in allGen if
            gen.bus.number == bus && occursin("Trip", gen.name) == false
        ])
    end

    set_units_base_system!(sys, "DEVICE_BASE")
    generators = [g for g in get_components(Generator, sys)]

    replace_gens = [g for g in generators if g.bus.number in ibrBus]
    for g in replace_gens
        if occursin("Trip", g.name) == false
            gen = get_component(ThermalStandard, sys, g.name)
            set_base_power!(gen, busCap[g.bus.number] * (1 - GF - Gf))
            set_base_power!(gen.dynamic_injector, busCap[g.bus.number] * (1 - GF - Gf))

            gen = get_component(GenericBattery, sys, join(["GF_Battery-", g.bus.number]))
            set_base_power!(gen, busCap[g.bus.number] * GF)
            set_base_power!(gen.dynamic_injector, busCap[g.bus.number] * GF)

            gen = get_component(GenericBattery, sys, join(["Gf_Battery-", g.bus.number]))
            set_base_power!(gen, busCap[g.bus.number] * Gf)
            set_base_power!(gen.dynamic_injector, busCap[g.bus.number] * Gf)
        end
    end
end

#= function change_ibr_penetration!(sys, GF, Gf, ibr_bus, bus_capacity)
    
    set_units_base_system!(sys, "DEVICE_BASE")
    generators = [g for g in get_components(Generator, sys)]
    
    replace_gens = [g for g in generators if g.bus.number in ibr_bus]
    
    for g in replace_gens
        if occursin("Trip", g.name)==false
            gen = get_component(ThermalStandard, sys, g.name)
            set_base_power!(gen, bus_capacity[g.bus.number]*(1-GF-Gf))
            set_base_power!(gen.dynamic_injector, bus_capacity[g.bus.number]*(1-GF-Gf))

            gen = get_component(GenericBattery, sys, join(["GF_Battery-", g.bus.number]))
            set_base_power!(gen, bus_capacity[g.bus.number]*GF)
            set_base_power!(gen.dynamic_injector, bus_capacity[g.bus.number]*GF)
                
            gen = get_component(GenericBattery, sys, join(["Gf_Battery-", g.bus.number]))
            set_base_power!(gen, bus_capacity[g.bus.number]*Gf)
            set_base_power!(gen.dynamic_injector, bus_capacity[g.bus.number]*Gf)
        end
    end

end

function getSystemProperties(sys)

    syncGen = collect(get_components(Generator, sys));
    ibrGen =  collect(get_components(GenericBattery, sys));
    allGen = vcat(syncGen, ibrGen);

    ibrBus = unique([gen.bus.number for gen in ibrGen]);
    syncBus = unique([gen.bus.number for gen in syncGen]);
    genBus = unique(vcat(ibrBus, syncBus));

    busCap=Dict(zip(genBus, zeros(length(genBus))))

    for bus in genBus
        busCap[bus] = sum([get_base_power(gen) for gen in allGen if gen.bus.number == bus && occursin("Trip", gen.name)==false ])
    end

    powerfow=solve_powerflow(sys)
    totalGen=sum(powerfow["bus_results"].P_gen)

    return busCap, totalGen, ibrBus, ibrGen, syncGen
end =#

#TODO - add a GFL 

#TODO - add a zoom option for visualize dataset (to see the fast transients)
#TODO - add P and Q to the plots
function visualize_dataset(dataset::Vector{PSIDS.SteadyStateNODEData})
    p1 = plot()
    p2 = plot()
    p3 = plot()
    p4 = plot()
    p5 = plot()
    p6 = plot()
    p7 = plot()
    p8 = plot()
    total_faults = length(dataset)
    for (ix, d) in enumerate(dataset)
        #@assert size(dataset[1].groundtruth_current)[1] == 2    #Surrogate should have a single connecting line 
        if dataset[ix].stable == true
            Vr = dataset[ix].surrogate_real_voltage[:]
            Vi = dataset[ix].surrogate_imag_voltage[:]
            Ir = dataset[ix].real_current[:]
            Ii = dataset[ix].imag_current[:]
            Vm = sqrt.(Vr .^ 2 + Vi .^ 2)
            Im = sqrt.(Ir .^ 2 + Ii .^ 2)
            Vθ = atan.(Vi ./ Vr)
            Iθ = atan.(Ii ./ Ir)
            plot!(p1, dataset[ix].tsteps, Vr, title = "Vr", legend = false)
            plot!(p2, dataset[ix].tsteps, Ir, title = "Ir", legend = false)
            plot!(p3, dataset[ix].tsteps, Vi, title = "Vi", legend = false)
            plot!(p4, dataset[ix].tsteps, Ii, title = "Ii", legend = false)
            plot!(p5, dataset[ix].tsteps, Vm, title = "Vm", legend = false)
            plot!(p6, dataset[ix].tsteps, Im, title = "Im", legend = false)
            plot!(p7, dataset[ix].tsteps, Vθ, title = "Vθ", legend = false)
            plot!(p8, dataset[ix].tsteps, Iθ, title = "Iθ", legend = false)
        end
    end
    plot!(p8, label = "test")
    total_faults = length(dataset)
    stable_faults = length(filter(x -> x.stable == true, dataset))
    unstable_faults = length(filter(x -> x.stable == false, dataset))
    @warn "\n connecting resistance: $(dataset[1].connecting_resistance)\n"
    @warn "\n connecting reactance: $(dataset[1].connecting_reactance)\n"
    @warn "\n total faults: $total_faults \n stable faults: $stable_faults \n unstable faults: $unstable_faults \n"
    l_tstops = [length(d.tstops) for d in dataset]
    l_tsteps = [length(d.tsteps) for d in dataset]

    @warn "length tstops", l_tstops
    @warn "length tsteps", l_tsteps

    return plot(p1, p2, p3, p4, p5, p6, p7, p8, layout = (4, 2), size = (700, 800))
end

function add_battery(sys, battery_name, bus_name, capacity, P, Q)
    return GenericBattery(
        name = battery_name,
        bus = get_component(Bus, sys, bus_name),
        available = true,
        prime_mover = PrimeMovers.BA,
        active_power = P,
        reactive_power = Q,
        rating = 1.1,
        base_power = capacity,
        initial_energy = 50.0,
        state_of_charge_limits = (min = 5.0, max = 100.0),
        input_active_power_limits = (min = 0.0, max = 1.0),
        output_active_power_limits = (min = 0.0, max = 1.0),
        reactive_power_limits = (min = -1.0, max = 1.0),
        efficiency = (in = 0.80, out = 0.90),
    )
end
