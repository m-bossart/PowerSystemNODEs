include("InverterModels.jl")
function get_init_gfm(p, vr_filter, vi_filter, ir_filter, ii_filter)
    return  (dx, x) -> begin
            dx[1] = 0
            x[1] = vi_filter
            dx[10]
            x[10] = vr_filter
            dx[9] = 0
            x[9] = ir_filter
            dx[19]
            x[19] = ii_filter
            gfm(dx,x,p,0)
            nothing
        end
end

function get_init_gfm_mm(p, vr_filter, vi_filter, ir_filter, ii_filter)
    return (dx, x) -> begin
                dx[1] = 0
                x[1] = vi_filter
                dx[10]
                x[10] = vr_filter
                dx[9] = 0
                x[9] = ir_filter
                dx[19]
                x[19] = ii_filter
                gfm_nn_out(dx,x,p,0)
                nothing
        end
end
