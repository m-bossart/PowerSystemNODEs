include("InverterModels.jl")
function get_init_gfm(p, ir_filter, ii_filter)
    return  (dx, x) -> begin
            dx[5] = 0
            x[5] = ir_filter
            dx[19] = 0
            x[19] = ii_filter
            gfm(dx,x,p,0)
            nothing
        end
end

#This fixes states 5 and 19, but in reality I want to fix
#the sum of 5 and 19 to the given value...
#Need to add two alebraic states which are the
#total real and reactive currents out.
#Change the mass matrix, etc....

function get_init_gfm_nn_states(p, ir_filter, ii_filter)
    return (dx, x) -> begin
                dx[24] = 0
                x[24] = ir_filter
                dx[25] = 0
                x[25] = ii_filter
                gfm_nn_states(dx,x,p,0)
                nothing
        end
end

function get_init_gfm_nn_voltage(p, ir_filter, ii_filter)
    return (dx, x) -> begin
                dx[24] = 0
                x[24] = ir_filter
                dx[25] = 0
                x[25] = ii_filter
                gfm_nn_voltage(dx,x,p,0)
                nothing
        end
end
