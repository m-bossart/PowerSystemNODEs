include("InverterModels.jl")




function get_init_gfm(p, ir_filter, ii_filter)
    return  (dx, x) -> begin
            dx[5] = 0
            x[5] = ir_filter
            dx[19] = 0
            x[19] = ii_filter
            gfm(dx,x,p,0.0)
            nothing
        end
end

function get_init_gfm_nn(p, ir_filter, ii_filter)
    return (dx, x) -> begin
                dx[24] = 0
                x[24] = ir_filter
                dx[25] = 0
                x[25] = ii_filter
                gfm_nn(dx,x,p,0.0)
                nothing
        end
end


function get_init_gfm_nn_states(p, ir_filter, ii_filter)
    return (dx, x) -> begin
                dx[24] = 0
                x[24] = ir_filter
                dx[25] = 0
                x[25] = ii_filter
                gfm_nn_states(dx,x,p,0.0)
                nothing
        end
end
