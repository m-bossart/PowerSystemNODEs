using InfrastructureSystems
const IS = InfrastructureSystems
using PowerSystems

s1 = System(100.0)
add_component!(s1,Bus(nothing))
source = Source(nothing)
pvs = PeriodicVariableSource(nothing)
add_component!(s1, source)
add_component!(s1, pvs, source)
@assert length(collect(get_components(PeriodicVariableSource,s1))) !== 0

found_source = collect(get_components(Source, s1))[1]
found_pvs = get_dynamic_injector(found_source)
remove_component!(s1,found_pvs)
remove_component!(s1,found_source)

s2 = System(100.0)
add_component!(s2,Bus(nothing))
add_component!(s2, found_source)
add_component!(s2, found_pvs, found_source)
@assert length(collect(get_components(PeriodicVariableSource,s2))) !== 0
