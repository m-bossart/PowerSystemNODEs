#Naming of surrogate models is as follows: <analytical_model>_<nn_external_inputs>_<feed_back_current>_<#_feedback_states>
#<analytical_model>: the odes
# vsm (virtual synchrnous machine)
#<nn_external_inputs>: the inputs to the model that are not states of the nn
# v: terminal voltage (real and imag)
# i: output current of the analytical model (real and imag)
# NOTE: could explore using states of the analytical model in the future
#<#_feedback_states>: number of outputs of the nn. The states associated with each output are automatically inputs to the nn
# 1-10
# NOTE: Input dimension of nn =  <nn_external_inputs> + <#_feedback_states>

mutable struct SurrParams
    nn::Vector{}                    #needs to accept dual numbers during training 
    ode::Vector{Float64}
    pf::Vector{Float64}
    network::Vector{Float64}
    scale::Vector{Float64}
    n_weights::Vector{Float64}
end

function vectorize(P::SurrParams)
    return vcat(P.nn, P.ode, P.pf, P.network, P.scale, P.n_weights)
end

function SurrParams(;)
    SurrParams([], [], [], [], [], [])
end

function none_v_t_0(dx, x, p, t, nn, Vm, Vθ)
    #PARAMETERS
    n_weights_nn = p[end]
    p_nn = p[Int(1):Int(n_weights_nn)]
    p_fixed = p[(Int(n_weights_nn) + 1):(end - 1)]

    P_pf = p_fixed[1]
    Q_pf = p_fixed[2]
    V_pf = p_fixed[3]
    θ_pf = p_fixed[4]
    Xtrans = p_fixed[5]
    Rtrans = p_fixed[6]
    V_scale = p_fixed[7]
    nn_scale = p_fixed[8]

    #STATE INDEX AND STATES
    i__ir_nn, ir_nn = 1, x[1]
    i__ii_nn, ii_nn = 2, x[2]

    Vr_pcc = Vm(t) * cos(Vθ(t)) + (ir_nn * Rtrans - ii_nn * Xtrans)
    Vi_pcc = Vm(t) * sin(Vθ(t)) + (ir_nn * Xtrans + ii_nn * Rtrans)

    nn_input = [
        P_pf,
        Q_pf,
        V_pf,
        θ_pf,
        (Vr_pcc - V_pf * cos(θ_pf)) * V_scale,
        (Vi_pcc - V_pf * sin(θ_pf)) * V_scale,
        ir_nn,
        ii_nn,
    ]

    #NODE   
    dx[i__ir_nn] = nn(nn_input, p_nn)[1] * nn_scale
    dx[i__ii_nn] = nn(nn_input, p_nn)[2] * nn_scale
end

function none_v_t_1(dx, x, p, t, nn, Vm, Vθ)
    #PARAMETERS
    n_weights_nn = p[end]
    p_nn = p[Int(1):Int(n_weights_nn)]
    p_fixed = p[(Int(n_weights_nn) + 1):(end - 1)]

    P_pf = p_fixed[1]
    Q_pf = p_fixed[2]
    V_pf = p_fixed[3]
    θ_pf = p_fixed[4]
    Xtrans = p_fixed[5]
    Rtrans = p_fixed[6]
    V_scale = p_fixed[7]
    nn_scale = p_fixed[8]

    #STATE INDEX AND STATES

    i__ir_nn, ir_nn = 1, x[1]
    i__ii_nn, ii_nn = 2, x[2]
    i__fb1, fb1 = 3, x[3]

    Vr_pcc = Vm(t) * cos(Vθ(t)) + (ir_nn * Rtrans - ii_nn * Xtrans)
    Vi_pcc = Vm(t) * sin(Vθ(t)) + (ir_nn * Xtrans + ii_nn * Rtrans)

    nn_input = [
        P_pf,
        Q_pf,
        V_pf,
        θ_pf,
        (Vr_pcc - V_pf * cos(θ_pf)) * V_scale,
        (Vi_pcc - V_pf * sin(θ_pf)) * V_scale,
        ir_nn,
        ii_nn,
        fb1,
    ]

    #NODE   
    dx[i__ir_nn] = nn(nn_input, p_nn)[1] * nn_scale
    dx[i__ii_nn] = nn(nn_input, p_nn)[2] * nn_scale
    dx[i__fb1] = nn(nn_input, p_nn)[3]
end

function none_v_t_2(dx, x, p, t, nn, Vm, Vθ)
    #PARAMETERS
    n_weights_nn = p[end]
    p_nn = p[Int(1):Int(n_weights_nn)]
    p_fixed = p[(Int(n_weights_nn) + 1):(end - 1)]

    P_pf = p_fixed[1]
    Q_pf = p_fixed[2]
    V_pf = p_fixed[3]
    θ_pf = p_fixed[4]
    Xtrans = p_fixed[5]
    Rtrans = p_fixed[6]
    V_scale = p_fixed[7]
    nn_scale = p_fixed[8]

    #STATE INDEX AND STATES

    i__ir_nn, ir_nn = 1, x[1]
    i__ii_nn, ii_nn = 2, x[2]
    i__fb1, fb1 = 3, x[3]
    i__fb2, fb2 = 4, x[4]

    Vr_pcc = Vm(t) * cos(Vθ(t)) + (ir_nn * Rtrans - ii_nn * Xtrans)
    Vi_pcc = Vm(t) * sin(Vθ(t)) + (ir_nn * Xtrans + ii_nn * Rtrans)

    nn_input = [
        P_pf,
        Q_pf,
        V_pf,
        θ_pf,
        (Vr_pcc - V_pf * cos(θ_pf)) * V_scale,
        (Vi_pcc - V_pf * sin(θ_pf)) * V_scale,
        ir_nn,
        ii_nn,
        fb1,
        fb2,
    ]

    #NODE   
    dx[i__ir_nn] = nn(nn_input, p_nn)[1] * nn_scale
    dx[i__ii_nn] = nn(nn_input, p_nn)[2] * nn_scale
    dx[i__fb1] = nn(nn_input, p_nn)[3]
    dx[i__fb2] = nn(nn_input, p_nn)[4]
end

function none_v_t_3(dx, x, p, t, nn, Vm, Vθ)
    #PARAMETERS
    n_weights_nn = p[end]
    p_nn = p[Int(1):Int(n_weights_nn)]
    p_fixed = p[(Int(n_weights_nn) + 1):(end - 1)]

    P_pf = p_fixed[1]
    Q_pf = p_fixed[2]
    V_pf = p_fixed[3]
    θ_pf = p_fixed[4]
    Xtrans = p_fixed[5]
    Rtrans = p_fixed[6]
    V_scale = p_fixed[7]
    nn_scale = p_fixed[8]

    #STATE INDEX AND STATES

    i__ir_nn, ir_nn = 1, x[1]
    i__ii_nn, ii_nn = 2, x[2]
    i__fb1, fb1 = 3, x[3]
    i__fb2, fb2 = 4, x[4]
    i__fb3, fb3 = 5, x[5]

    Vr_pcc = Vm(t) * cos(Vθ(t)) + (ir_nn * Rtrans - ii_nn * Xtrans)
    Vi_pcc = Vm(t) * sin(Vθ(t)) + (ir_nn * Xtrans + ii_nn * Rtrans)

    nn_input = [
        P_pf,
        Q_pf,
        V_pf,
        θ_pf,
        (Vr_pcc - V_pf * cos(θ_pf)) * V_scale,
        (Vi_pcc - V_pf * sin(θ_pf)) * V_scale,
        ir_nn,
        ii_nn,
        fb1,
        fb2,
        fb3,
    ]

    #NODE   
    dx[i__ir_nn] = nn(nn_input, p_nn)[1] * nn_scale
    dx[i__ii_nn] = nn(nn_input, p_nn)[2] * nn_scale
    dx[i__fb1] = nn(nn_input, p_nn)[3]
    dx[i__fb2] = nn(nn_input, p_nn)[4]
    dx[i__fb3] = nn(nn_input, p_nn)[5]
end

function none_v_t_4(dx, x, p, t, nn, Vm, Vθ)
    #PARAMETERS
    n_weights_nn = p[end]
    p_nn = p[Int(1):Int(n_weights_nn)]
    p_fixed = p[(Int(n_weights_nn) + 1):(end - 1)]

    P_pf = p_fixed[1]
    Q_pf = p_fixed[2]
    V_pf = p_fixed[3]
    θ_pf = p_fixed[4]
    Xtrans = p_fixed[5]
    Rtrans = p_fixed[6]
    V_scale = p_fixed[7]
    nn_scale = p_fixed[8]

    #STATE INDEX AND STATES

    i__ir_nn, ir_nn = 1, x[1]
    i__ii_nn, ii_nn = 2, x[2]
    i__fb1, fb1 = 3, x[3]
    i__fb2, fb2 = 4, x[4]
    i__fb3, fb3 = 5, x[5]
    i__fb4, fb4 = 6, x[6]

    Vr_pcc = Vm(t) * cos(Vθ(t)) + (ir_nn * Rtrans - ii_nn * Xtrans)
    Vi_pcc = Vm(t) * sin(Vθ(t)) + (ir_nn * Xtrans + ii_nn * Rtrans)

    nn_input = [
        P_pf,
        Q_pf,
        V_pf,
        θ_pf,
        (Vr_pcc - V_pf * cos(θ_pf)) * V_scale,
        (Vi_pcc - V_pf * sin(θ_pf)) * V_scale,
        ir_nn,
        ii_nn,
        fb1,
        fb2,
        fb3,
        fb4,
    ]

    #NODE   
    dx[i__ir_nn] = nn(nn_input, p_nn)[1] * nn_scale
    dx[i__ii_nn] = nn(nn_input, p_nn)[2] * nn_scale
    dx[i__fb1] = nn(nn_input, p_nn)[3]
    dx[i__fb2] = nn(nn_input, p_nn)[4]
    dx[i__fb3] = nn(nn_input, p_nn)[5]
    dx[i__fb4] = nn(nn_input, p_nn)[6]
end

function none_v_t_5(dx, x, p, t, nn, Vm, Vθ)
    #PARAMETERS
    n_weights_nn = p[end]
    p_nn = p[Int(1):Int(n_weights_nn)]
    p_fixed = p[(Int(n_weights_nn) + 1):(end - 1)]

    P_pf = p_fixed[1]
    Q_pf = p_fixed[2]
    V_pf = p_fixed[3]
    θ_pf = p_fixed[4]
    Xtrans = p_fixed[5]
    Rtrans = p_fixed[6]
    V_scale = p_fixed[7]
    nn_scale = p_fixed[8]

    #STATE INDEX AND STATES

    i__ir_nn, ir_nn = 1, x[1]
    i__ii_nn, ii_nn = 2, x[2]
    i__fb1, fb1 = 3, x[3]
    i__fb2, fb2 = 4, x[4]
    i__fb3, fb3 = 5, x[5]
    i__fb4, fb4 = 6, x[6]
    i__fb5, fb5 = 7, x[7]

    Vr_pcc = Vm(t) * cos(Vθ(t)) + (ir_nn * Rtrans - ii_nn * Xtrans)
    Vi_pcc = Vm(t) * sin(Vθ(t)) + (ir_nn * Xtrans + ii_nn * Rtrans)

    nn_input = [
        P_pf,
        Q_pf,
        V_pf,
        θ_pf,
        (Vr_pcc - V_pf * cos(θ_pf)) * V_scale,
        (Vi_pcc - V_pf * sin(θ_pf)) * V_scale,
        ir_nn,
        ii_nn,
        fb1,
        fb2,
        fb3,
        fb4,
        fb5,
    ]

    #NODE   
    dx[i__ir_nn] = nn(nn_input, p_nn)[1] * nn_scale
    dx[i__ii_nn] = nn(nn_input, p_nn)[2] * nn_scale
    dx[i__fb1] = nn(nn_input, p_nn)[3]
    dx[i__fb2] = nn(nn_input, p_nn)[4]
    dx[i__fb3] = nn(nn_input, p_nn)[5]
    dx[i__fb4] = nn(nn_input, p_nn)[6]
    dx[i__fb5] = nn(nn_input, p_nn)[7]
end

function vsm_v_t_0(dx, x, p, t, nn, Vm, Vθ)
    #PARAMETERS
    n_weights_nn = p[end]
    p_nn = p[Int(1):Int(n_weights_nn)]
    p_fixed = p[(Int(n_weights_nn) + 1):(end - 1)]

    ω_lp = p_fixed[1]
    kp_pll = p_fixed[2]
    ki_pll = p_fixed[3]
    Ta = p_fixed[4]
    kd = p_fixed[5]
    kω = p_fixed[6]
    kq = p_fixed[7]
    ωf = p_fixed[8]
    kpv = p_fixed[9]
    kiv = p_fixed[10]
    kffv = p_fixed[11]
    rv = p_fixed[12]
    lv = p_fixed[13]
    kpc = p_fixed[14]
    kic = p_fixed[15]
    kffi = p_fixed[16]
    ωad = p_fixed[17]
    kad = p_fixed[18]
    lf = p_fixed[19]
    rf = p_fixed[20]
    cf = p_fixed[21]
    lg = p_fixed[22]
    rg = p_fixed[23]
    Vref = p_fixed[24]
    ωref = p_fixed[25]
    Pref = p_fixed[26]
    Qref = p_fixed[27]
    P_pf = p_fixed[28]
    Q_pf = p_fixed[29]
    V_pf = p_fixed[30]
    θ_pf = p_fixed[31]
    Xtrans = p_fixed[32]
    Rtrans = p_fixed[33]
    V_scale = p_fixed[34]
    nn_scale = p_fixed[35]

    #STATE INDEX AND STATES
    i__ir_out, ir_out = 1, x[1]
    i__ii_out, ii_out = 2, x[2]
    i__ir_nn, ir_nn = 3, x[3]
    i__ii_nn, ii_nn = 4, x[4]
    i__vi_filter, vi_filter = 5, x[5]
    i__γd_ic, γd_ic = 6, x[6]
    i__vq_pll, vq_pll = 7, x[7]
    i__γq_ic, γq_ic = 8, x[8]
    i__ir_filter, ir_filter = 9, x[9]
    i__ξd_ic, ξd_ic = 10, x[10]
    i__ϕd_ic, ϕd_ic = 11, x[11]
    i__ε_pll, ε_pll = 12, x[12]
    i__ir_cnv, ir_cnv = 13, x[13]
    i__vr_filter, vr_filter = 14, x[14]
    i__ω_oc, ω_oc = 15, x[15]
    i__ξq_ic, ξq_ic = 16, x[16]
    i__vd_pll, vd_pll = 17, x[17]
    i__q_oc, q_oc = 18, x[18]
    i__ϕq_ic, ϕq_ic = 19, x[19]
    i__θ_pll, θ_pll = 20, x[20]
    i__θ_oc, θ_oc = 21, x[21]
    i__ii_cnv, ii_cnv = 22, x[22]
    i__ii_filter, ii_filter = 23, x[23]

    ω_base = 60.0 * 2 * pi
    ω_sys = 1.0

    #PLL
    δω_pll = kp_pll * atan(vq_pll / vd_pll) + ki_pll * ε_pll
    ω_pll = δω_pll + ω_sys
    vd_filt_pll = sin(θ_pll + pi / 2) * vr_filter - cos(θ_pll + pi / 2) * vi_filter
    vq_filt_pll = cos(θ_pll + pi / 2) * vr_filter + sin(θ_pll + pi / 2) * vi_filter

    dx[i__vd_pll] = ω_lp * (vd_filt_pll - vd_pll)                                 #docs:(1a)
    dx[i__vq_pll] = ω_lp * (vq_filt_pll - vq_pll)                                 #docs:(1b)
    dx[i__ε_pll] = atan(vq_pll / vd_pll)                                         #docs:(1c)
    dx[i__θ_pll] = ω_base * δω_pll                                             #docs:(1d)

    #OUTER LOOP CONTROL
    pe = vr_filter * ir_filter + vi_filter * ii_filter                              #docs:(2d)
    qe = vi_filter * ir_filter - vr_filter * ii_filter                              #docs:(2e)
    v_ref_olc = Vref + kq * (Qref - q_oc)

    dx[i__ω_oc] = (Pref - pe - kd * (ω_oc - ω_pll) - kω * (ω_oc - ωref)) / Ta        #docs:(2a)
    dx[i__θ_oc] = ω_base * (ω_oc - ω_sys)                                            #docs:(2b)
    dx[i__q_oc] = ωf * (qe - q_oc)                                               #docs:(2c)

    #INNER LOOP CONTROL
    #reference transormations
    vd_filt_olc = sin(θ_oc + pi / 2) * vr_filter - cos(θ_oc + pi / 2) * vi_filter
    vq_filt_olc = cos(θ_oc + pi / 2) * vr_filter + sin(θ_oc + pi / 2) * vi_filter
    id_filt_olc = sin(θ_oc + pi / 2) * ir_filter - cos(θ_oc + pi / 2) * ii_filter
    iq_filt_olc = cos(θ_oc + pi / 2) * ir_filter + sin(θ_oc + pi / 2) * ii_filter
    id_cnv_olc = sin(θ_oc + pi / 2) * ir_cnv - cos(θ_oc + pi / 2) * ii_cnv
    iq_cnv_olc = cos(θ_oc + pi / 2) * ir_cnv + sin(θ_oc + pi / 2) * ii_cnv

    #Voltage control equations
    Vd_filter_ref = v_ref_olc - rv * id_filt_olc + ω_oc * lv * iq_filt_olc      #docs:(3g)
    Vq_filter_ref = -rv * iq_filt_olc - ω_oc * lv * id_filt_olc                 #docs:(3h)
    dx[i__ξd_ic] = Vd_filter_ref - vd_filt_olc                                 #docs:(3a)
    dx[i__ξq_ic] = Vq_filter_ref - vq_filt_olc                                 #docs:(3b)

    #current control equations
    Id_cnv_ref = (
        kpv * (Vd_filter_ref - vd_filt_olc) + kiv * ξd_ic -                           #docs:(3i)
        cf * ω_oc * vq_filt_olc + kffi * id_filt_olc
    )
    Iq_cnv_ref = (
        kpv * (Vq_filter_ref - vq_filt_olc) +
        kiv * ξq_ic +                           #docs:(3j)
        cf * ω_oc * vd_filt_olc +
        kffi * iq_filt_olc
    )
    dx[i__γd_ic] = Id_cnv_ref - id_cnv_olc                                      #docs:(3c)
    dx[i__γq_ic] = Iq_cnv_ref - iq_cnv_olc                                      #docs:(3d)

    #active damping equations
    Vd_cnv_ref = (
        kpc * (Id_cnv_ref - id_cnv_olc) + kic * γd_ic - lf * ω_oc * iq_cnv_olc +          #docs:(3k)
        kffv * vd_filt_olc - kad * (vd_filt_olc - ϕd_ic)
    )
    Vq_cnv_ref = (
        kpc * (Iq_cnv_ref - iq_cnv_olc) +
        kic * γq_ic +
        lf * ω_oc * id_cnv_olc +          #docs:(3l)
        kffv * vq_filt_olc - kad * (vq_filt_olc - ϕq_ic)
    )
    dx[i__ϕd_ic] = ωad * (vd_filt_olc - ϕd_ic)                                 #docs:(3e)
    dx[i__ϕq_ic] = ωad * (vq_filt_olc - ϕq_ic)                                  #docs:(3f)

    #LCL FILTER
    #reference transformations
    Vr_cnv = sin(θ_oc + pi / 2) * Vd_cnv_ref + cos(θ_oc + pi / 2) * Vq_cnv_ref
    Vi_cnv = -cos(θ_oc + pi / 2) * Vd_cnv_ref + sin(θ_oc + pi / 2) * Vq_cnv_ref

    Vr_pcc = Vm(t) * cos(Vθ(t)) + (ir_out * Rtrans - ii_out * Xtrans)
    Vi_pcc = Vm(t) * sin(Vθ(t)) + (ir_out * Xtrans + ii_out * Rtrans)

    dx[i__ir_cnv] =
        (ω_base / lf) *                                                #docs:(5a)
        (Vr_cnv - vr_filter - rf * ir_cnv + ω_sys * lf * ii_cnv)
    dx[i__ii_cnv] =
        (ω_base / lf) *                                                #docs:(5b)
        (Vi_cnv - vi_filter - rf * ii_cnv - ω_sys * lf * ir_cnv)
    dx[i__vr_filter] =
        (ω_base / cf) *                                             #docs:(5c)
        (ir_cnv - ir_filter + ω_sys * cf * vi_filter)
    dx[i__vi_filter] =
        (ω_base / cf) *                                             #docs:(5d)
        (ii_cnv - ii_filter - ω_sys * cf * vr_filter)
    dx[i__ir_filter] =
        (ω_base / lg) *                                             #docs:(5e)
        (vr_filter - Vr_pcc - rg * ir_filter + ω_sys * lg * ii_filter)
    dx[i__ii_filter] =
        +(ω_base / lg) *                                            #docs:(5f)
        (vi_filter - Vi_pcc - rg * ii_filter - ω_sys * lg * ir_filter)
    #NN INPUT
    Vr0 = V_pf * cos(θ_pf)
    Vi0 = V_pf * sin(θ_pf)
    nn_input = [
        P_pf,
        Q_pf,
        V_pf,
        θ_pf,
        (Vr_pcc - Vr0) * V_scale,
        (Vi_pcc - Vi0) * V_scale,
        ir_nn,
        ii_nn,
    ]

    #NN CURRENT SOURCE
    dx[i__ir_nn] = nn(nn_input, p_nn)[1] * nn_scale
    dx[i__ii_nn] = nn(nn_input, p_nn)[2] * nn_scale

    #ALEGRAIC STATE - OUTPUT CURRENT
    dx[i__ir_out] = ir_out - ir_nn - ir_filter
    dx[i__ii_out] = ii_out - ii_nn - ii_filter
end
