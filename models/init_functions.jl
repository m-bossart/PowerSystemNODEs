include("SurrogateModels.jl")




function get_init_vsm_nn_v_2(p, ir_filter, ii_filter)
    return (dx, x) -> begin
                dx[22] = 0
                x[22] = ir_filter
                dx[23] = 0
                x[23] = ii_filter
                vsm_nn_v_2(dx,x,p,0.0)
                nothing
        end
end

function get_init_vsm_nn_v_3(p, ir_filter, ii_filter)
        return (dx, x) -> begin
                    dx[22] = 0
                    x[22] = ir_filter
                    dx[23] = 0
                    x[23] = ii_filter
                    vsm_nn_v_3(dx,x,p,0.0)
                    nothing
        end
end
    

function get_init_vsm_nn_v_4(p, ir_filter, ii_filter)
        return (dx, x) -> begin
                    dx[22] = 0
                    x[22] = ir_filter
                    dx[23] = 0
                    x[23] = ii_filter
                    vsm_nn_v_4(dx,x,p,0.0)
                    nothing
        end
end

function get_init_vsm_nn_v_5(p, ir_filter, ii_filter)
        return (dx, x) -> begin
                    dx[22] = 0
                    x[22] = ir_filter
                    dx[23] = 0
                    x[23] = ii_filter
                    vsm_nn_v_5(dx,x,p,0.0)
                    nothing
        end
end
    


    

    
