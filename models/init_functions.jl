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
#THURSDAY TO DO!
function get_init_gfm_nn(p, ir_filter, ii_filter)
    return (dx, x) -> begin
                dx[5] = 0
                x[5] = ir_filter
                dx[19] = 0
                x[19] = ii_filter
                gfm_nn(dx,x,p,0)
                nothing
        end
end
