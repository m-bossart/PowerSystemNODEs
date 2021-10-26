include("SurrogateModels.jl")

function get_init_vsm_v_t_0(p, ir_filter, ii_filter, nn, Vm, Vθ)
    return (dx, x) -> begin
                dx[22] = 0
                x[22] = ir_filter
                dx[23] = 0
                x[23] = ii_filter
                vsm_v_t_0(dx,x,p,0.0, nn, Vm, Vθ)
                nothing
        end
end
