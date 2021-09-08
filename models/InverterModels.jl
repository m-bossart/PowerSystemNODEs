function gfm(dx,x,p,t)
    #PARAMETERS
    ω_lp = p[1]
    kp_pll = p[2]
    ki_pll = p[3]
    Ta = p[4]
    kd = p[5]
    kω = p[6]
    kq = p[7]
    ωf = p[8]
    kpv = p[9]
    kiv = p[10]
    kffv = p[11]
    rv = p[12]
    lv = p[13]
    kpc = p[14]
    kic = p[15]
    kffi = p[16]
    ωad = p[17]
    kad = p[18]
    lf = p[19]
    rf = p[20]
    cf = p[21]
    lg = p[22]
    rg = p[23]
    Vref = p[24]    #Reference treated as parameter, but should NOT get updated during training
    ωref = p[25]    #Reference treated as parameter, but should NOT get updated during training
    Pref = p[26]    #Reference treated as parameter, but should NOT get updated during training
    Qref = p[27]    #Reference treated as parameter, but should NOT get updated during training
    Xtrans = p[28]
    Rtrans = p[29]

    #STATE INDEX AND STATES
    i__vi_filter, vi_filter = 1, x[1]
    i__γd_ic, γd_ic = 2, x[2]
    i__vq_pll, vq_pll = 3, x[3]
    i__γq_ic, γq_ic = 4, x[4]
    i__ir_filter, ir_filter = 5, x[5]
    i__ξd_ic, ξd_ic = 6, x[6]
    i__ϕd_ic, ϕd_ic = 7, x[7]
    i__ε_pll, ε_pll = 8, x[8]
    i__ir_cnv, ir_cnv = 9, x[9]
    i__vr_filter, vr_filter = 10, x[10]
    i__ω_oc, ω_oc = 11, x[11]
    i__ξq_ic, ξq_ic = 12, x[12]
    i__vd_pll, vd_pll = 13, x[13]
    i__q_oc, q_oc = 14, x[14]
    i__ϕq_ic, ϕq_ic = 15, x[15]
    i__θ_pll, θ_pll = 16, x[16]
    i__θ_oc, θ_oc = 17, x[17]
    i__ii_cnv, ii_cnv = 18, x[18]
    i__ii_filter, ii_filter = 19, x[19]

    ω_base = 60.0*2*pi
    ω_sys = 1.0

    #PLL
    δω_pll = kp_pll * atan(vq_pll/vd_pll) + ki_pll* ε_pll
    ω_pll = δω_pll + ω_sys
    vd_filt_pll = sin(θ_pll + pi/2)*vr_filter - cos(θ_pll + pi/2)*vi_filter
    vq_filt_pll = cos(θ_pll + pi/2)*vr_filter + sin(θ_pll + pi/2)*vi_filter

    dx[i__vd_pll]=  ω_lp*(vd_filt_pll - vd_pll)              #docs:(1a)
    dx[i__vq_pll]=  ω_lp*(vq_filt_pll - vq_pll)              #docs:(1b)
    dx[i__ε_pll] =  atan(vq_pll/vd_pll)                      #docs:(1c)
    dx[i__θ_pll] =  ω_base * δω_pll                          #docs:(1d)

    #OUTER LOOP CONTROL
    pe = vr_filter*ir_filter + vi_filter*ii_filter           #docs:(2d)
    qe = vi_filter*ir_filter - vr_filter*ii_filter           #docs:(2e)
    v_ref_olc = Vref + kq * (Qref - q_oc)


    dx[i__ω_oc]= (Pref - pe - kd*(ω_oc - ω_pll) - kω*(ω_oc - ωref)) / Ta #docs:(2a)
    dx[i__θ_oc]= ω_base*(ω_oc-ω_sys)                         #docs:(2b)
    dx[i__q_oc]= ωf * (qe - q_oc)                            #docs:(2c)

    #INNER LOOP CONTROL
    #reference transormations
    vd_filt_olc = sin(θ_oc + pi/2)*vr_filter - cos(θ_oc + pi/2)*vi_filter
    vq_filt_olc = cos(θ_oc + pi/2)*vr_filter + sin(θ_oc + pi/2)*vi_filter
    id_filt_olc = sin(θ_oc + pi/2)*ir_filter - cos(θ_oc + pi/2)*ii_filter
    iq_filt_olc = cos(θ_oc + pi/2)*ir_filter + sin(θ_oc + pi/2)*ii_filter
    id_cnv_olc =  sin(θ_oc + pi/2)*ir_cnv - cos(θ_oc + pi/2)*ii_cnv
    iq_cnv_olc =  cos(θ_oc + pi/2)*ir_cnv + sin(θ_oc + pi/2)*ii_cnv

    #Voltage control equations
    Vd_filter_ref = v_ref_olc - rv * id_filt_olc + ω_oc * lv * iq_filt_olc    #docs:(3g)
    Vq_filter_ref = -rv * iq_filt_olc - ω_oc * lv * id_filt_olc               #docs:(3h)
    dx[i__ξd_ic]=  Vd_filter_ref -  vd_filt_olc               #docs:(3a)
    dx[i__ξq_ic]=  Vq_filter_ref -  vq_filt_olc               #docs:(3b)

    #current control equations
    Id_cnv_ref =(
      kpv*(Vd_filter_ref - vd_filt_olc) + kiv*ξd_ic -         #docs:(3i)
      cf*ω_oc*vq_filt_olc + kffi*id_filt_olc
    )
    Iq_cnv_ref =(
      kpv*(Vq_filter_ref - vq_filt_olc) + kiv*ξq_ic +         #docs:(3j)
      cf*ω_oc*vd_filt_olc + kffi*iq_filt_olc
    )
    dx[i__γd_ic]=  Id_cnv_ref - id_cnv_olc                    #docs:(3c)
    dx[i__γq_ic]=  Iq_cnv_ref - iq_cnv_olc                    #docs:(3d)

    #active damping equations
    Vd_cnv_ref =(
      kpc*(Id_cnv_ref - id_cnv_olc) + kic*γd_ic - lf*ω_oc*iq_cnv_olc +        #docs:(3k)
      kffv*vd_filt_olc - kad*(vd_filt_olc-ϕd_ic)
    )
    Vq_cnv_ref =(
      kpc*(Iq_cnv_ref - iq_cnv_olc) + kic*γq_ic + lf*ω_oc*id_cnv_olc +        #docs:(3l)
      kffv*vq_filt_olc - kad*(vq_filt_olc-ϕq_ic)
    )
    dx[i__ϕd_ic] =  ωad * (vd_filt_olc - ϕd_ic)                #docs:(3e)
    dx[i__ϕq_ic] =  ωad * (vq_filt_olc - ϕq_ic)                #docs:(3f)

    #LCL FILTER
    #reference transformations
    Vr_cnv =   sin(θ_oc + pi/2)*Vd_cnv_ref + cos(θ_oc + pi/2)*Vq_cnv_ref
    Vi_cnv =  -cos(θ_oc + pi/2)*Vd_cnv_ref + sin(θ_oc + pi/2)*Vq_cnv_ref

    Vr_pcc = Vm(t)*cos(Vθ(t)) + (ir_filter*Rtrans - ii_filter*Xtrans)
    Vi_pcc = Vm(t)*sin(Vθ(t)) + (ir_filter*Xtrans + ii_filter*Rtrans)

    dx[i__ir_cnv]= (ω_base/lf) *                            #docs:(5a)
      (Vr_cnv - vr_filter - rf*ir_cnv + ω_sys*lf*ii_cnv)
    dx[i__ii_cnv]= (ω_base/lf) *                             #docs:(5b)
      (Vi_cnv - vi_filter - rf*ii_cnv - ω_sys*lf*ir_cnv)
    dx[i__vr_filter]= (ω_base/cf) *                       #docs:(5c)
      (ir_cnv - ir_filter + ω_sys*cf*vi_filter)
    dx[i__vi_filter]= (ω_base/cf) *                       #docs:(5d)
      (ii_cnv - ii_filter - ω_sys*cf*vr_filter)
    dx[i__ir_filter]= (ω_base/lg) *                       #docs:(5e)
      (vr_filter - Vr_pcc - rg*ir_filter + ω_sys*lg*ii_filter)
    dx[i__ii_filter]= +(ω_base/lg) *                       #docs:(5f)
      (vi_filter - Vi_pcc - rg*ii_filter - ω_sys*lg*ir_filter)
end

#19th order gfm model with nn current source on the output with input V(t)
#and output I(t)
function gfm_nn(dx,x,p,t)
    #PARAMETERS
    #n_weights_nn = Int(p[end])
    p_nn = p[1:n_weights_nn]
    p_ode = p[n_weights_nn+1:end]

    ω_lp = p_ode[1]
    kp_pll = p_ode[2]
    ki_pll = p_ode[3]
    Ta = p_ode[4]
    kd = p_ode[5]
    kω = p_ode[6]
    kq = p_ode[7]
    ωf = p_ode[8]
    kpv = p_ode[9]
    kiv = p_ode[10]
    kffv = p_ode[11]
    rv = p_ode[12]
    lv = p_ode[13]
    kpc = p_ode[14]
    kic = p_ode[15]
    kffi = p_ode[16]
    ωad = p_ode[17]
    kad = p_ode[18]
    lf = p_ode[19]
    rf = p_ode[20]
    cf = p_ode[21]
    lg = p_ode[22]
    rg = p_ode[23]
    Vref = p_ode[24]    #Reference treated as parameter, but should NOT get updated during training
    ωref = p_ode[25]    #Reference treated as parameter, but should NOT get updated during training
    Pref = p_ode[26]    #Reference treated as parameter, but should NOT get updated during training
    Qref = p_ode[27]    #Reference treated as parameter, but should NOT get updated during training
    Xtrans = p_ode[28]
    Rtrans = p_ode[29]
    ir_offset = p_ode[30]
    ii_offset = p_ode[31]

    #STATE INDEX AND STATES
    i__vi_filter, vi_filter = 1, x[1]
    i__γd_ic, γd_ic = 2, x[2]
    i__vq_pll, vq_pll = 3, x[3]
    i__γq_ic, γq_ic = 4, x[4]
    i__ir_filter, ir_filter = 5, x[5]
    i__ξd_ic, ξd_ic = 6, x[6]
    i__ϕd_ic, ϕd_ic = 7, x[7]
    i__ε_pll, ε_pll = 8, x[8]
    i__ir_cnv, ir_cnv = 9, x[9]
    i__vr_filter, vr_filter = 10, x[10]
    i__ω_oc, ω_oc = 11, x[11]
    i__ξq_ic, ξq_ic = 12, x[12]
    i__vd_pll, vd_pll = 13, x[13]
    i__q_oc, q_oc = 14, x[14]
    i__ϕq_ic, ϕq_ic = 15, x[15]
    i__θ_pll, θ_pll = 16, x[16]
    i__θ_oc, θ_oc = 17, x[17]
    i__ii_cnv, ii_cnv = 18, x[18]
    i__ii_filter, ii_filter = 19, x[19]
    i__ir_nn, ir_nn = 20, x[20]
    i__ii_nn, ii_nn = 21, x[21]
    #i__ir_offset, ir_offset = 22, x[22]
    #i__ii_offset, ii_offset =  23, x[23]
    i__ir_out, ir_out =  22, x[22]
    i__ii_out, ii_out =  23, x[23]

    ω_base = 60.0*2*pi
    ω_sys = 1.0

    #PLL
    δω_pll = kp_pll * atan(vq_pll/vd_pll) + ki_pll* ε_pll
    ω_pll = δω_pll + ω_sys
    vd_filt_pll = sin(θ_pll + pi/2)*vr_filter - cos(θ_pll + pi/2)*vi_filter
    vq_filt_pll = cos(θ_pll + pi/2)*vr_filter + sin(θ_pll + pi/2)*vi_filter

    dx[i__vd_pll]=  ω_lp*(vd_filt_pll - vd_pll)                                 #docs:(1a)
    dx[i__vq_pll]=  ω_lp*(vq_filt_pll - vq_pll)                                 #docs:(1b)
    dx[i__ε_pll] =  atan(vq_pll/vd_pll)                                         #docs:(1c)
    dx[i__θ_pll] =  ω_base * δω_pll                                             #docs:(1d)

    #OUTER LOOP CONTROL
    pe = vr_filter*ir_filter + vi_filter*ii_filter                              #docs:(2d)
    qe = vi_filter*ir_filter - vr_filter*ii_filter                              #docs:(2e)
    v_ref_olc = Vref + kq * (Qref - q_oc)


    dx[i__ω_oc]= (Pref - pe - kd*(ω_oc - ω_pll) - kω*(ω_oc - ωref)) / Ta        #docs:(2a)
    dx[i__θ_oc]= ω_base*(ω_oc-ω_sys)                                            #docs:(2b)
    dx[i__q_oc]= ωf * (qe - q_oc)                                               #docs:(2c)

    #INNER LOOP CONTROL
    #reference transormations
    vd_filt_olc = sin(θ_oc + pi/2)*vr_filter - cos(θ_oc + pi/2)*vi_filter
    vq_filt_olc = cos(θ_oc + pi/2)*vr_filter + sin(θ_oc + pi/2)*vi_filter
    id_filt_olc = sin(θ_oc + pi/2)*ir_filter - cos(θ_oc + pi/2)*ii_filter
    iq_filt_olc = cos(θ_oc + pi/2)*ir_filter + sin(θ_oc + pi/2)*ii_filter
    id_cnv_olc =  sin(θ_oc + pi/2)*ir_cnv - cos(θ_oc + pi/2)*ii_cnv
    iq_cnv_olc =  cos(θ_oc + pi/2)*ir_cnv + sin(θ_oc + pi/2)*ii_cnv

    #Voltage control equations
    Vd_filter_ref = v_ref_olc - rv * id_filt_olc + ω_oc * lv * iq_filt_olc      #docs:(3g)
    Vq_filter_ref = -rv * iq_filt_olc - ω_oc * lv * id_filt_olc                 #docs:(3h)
    dx[i__ξd_ic]=  Vd_filter_ref -  vd_filt_olc                                 #docs:(3a)
    dx[i__ξq_ic]=  Vq_filter_ref -  vq_filt_olc                                 #docs:(3b)

    #current control equations
    Id_cnv_ref =(
      kpv*(Vd_filter_ref - vd_filt_olc) + kiv*ξd_ic -                           #docs:(3i)
      cf*ω_oc*vq_filt_olc + kffi*id_filt_olc
    )
    Iq_cnv_ref =(
      kpv*(Vq_filter_ref - vq_filt_olc) + kiv*ξq_ic +                           #docs:(3j)
      cf*ω_oc*vd_filt_olc + kffi*iq_filt_olc
    )
    dx[i__γd_ic]=  Id_cnv_ref - id_cnv_olc                                      #docs:(3c)
    dx[i__γq_ic]=  Iq_cnv_ref - iq_cnv_olc                                      #docs:(3d)

    #active damping equations
    Vd_cnv_ref =(
      kpc*(Id_cnv_ref - id_cnv_olc) + kic*γd_ic - lf*ω_oc*iq_cnv_olc +          #docs:(3k)
      kffv*vd_filt_olc - kad*(vd_filt_olc-ϕd_ic)
    )
    Vq_cnv_ref =(
      kpc*(Iq_cnv_ref - iq_cnv_olc) + kic*γq_ic + lf*ω_oc*id_cnv_olc +          #docs:(3l)
      kffv*vq_filt_olc - kad*(vq_filt_olc-ϕq_ic)
    )
    dx[i__ϕd_ic] =  ωad * (vd_filt_olc - ϕd_ic)                                 #docs:(3e)
    dx[i__ϕq_ic]=  ωad * (vq_filt_olc - ϕq_ic)                                  #docs:(3f)

    #LCL FILTER
    #reference transformations
    Vr_cnv =   sin(θ_oc + pi/2)*Vd_cnv_ref + cos(θ_oc + pi/2)*Vq_cnv_ref
    Vi_cnv =  -cos(θ_oc + pi/2)*Vd_cnv_ref + sin(θ_oc + pi/2)*Vq_cnv_ref

    Vr_pcc = Vm(t)*cos(Vθ(t)) + (ir_out*Rtrans - ii_out*Xtrans)
    Vi_pcc = Vm(t)*sin(Vθ(t)) + (ir_out*Xtrans + ii_out*Rtrans)

    dx[i__ir_cnv]= (ω_base/lf) *                                                #docs:(5a)
      (Vr_cnv - vr_filter - rf*ir_cnv + ω_sys*lf*ii_cnv)
    dx[i__ii_cnv]= (ω_base/lf) *                                                #docs:(5b)
      (Vi_cnv - vi_filter - rf*ii_cnv - ω_sys*lf*ir_cnv)
    dx[i__vr_filter]= (ω_base/cf) *                                             #docs:(5c)
      (ir_cnv - ir_filter + ω_sys*cf*vi_filter)
    dx[i__vi_filter]= (ω_base/cf) *                                             #docs:(5d)
      (ii_cnv - ii_filter - ω_sys*cf*vr_filter)
    dx[i__ir_filter]= (ω_base/lg) *                                             #docs:(5e)
      (vr_filter - Vr_pcc - rg*ir_filter + ω_sys*lg*ii_filter)
    dx[i__ii_filter]= +(ω_base/lg) *                                            #docs:(5f)
      (vi_filter - Vi_pcc - rg*ii_filter - ω_sys*lg*ir_filter)

     #NN CURRENT SOURCE (NN input includes terminal voltage only)
     dx[i__ir_nn] =ir_offset -  nn([Vm(t), Vθ(t)], p_nn)[1] * nn_scale
     dx[i__ii_nn] =ii_offset -  nn([Vm(t), Vθ(t)], p_nn)[2] * nn_scale

     #Current offset for the NN Current Source
     #dx[i__ir_offset] = 0.0    #NOTE: Got an error using Zygote (reverse mode ad) if derivative is set to 0.0
     #dx[i__ii_offset] = 0.0    #NOTE: Got an error using Zygote (reverse mode ad) if derivative is set to 0.0

     #ALEGRAIC STATE - OUTPUT CURRENT
     dx[i__ir_out] =  ir_out - ir_nn - ir_filter
     dx[i__ii_out] =  ii_out - ii_nn - ii_filter
end



#19th order gfm model with nn current source on the output with input V(t)
#and output I(t)
"""
Surrogate with NN source depending on all 19 inverter states
"""
function gfm_nn_states(dx,x,p,t)
    #PARAMETERS
    p_nn = p[1:n_weights_nn_states]
    p_ode = p[n_weights_nn_states+1:end]

    ω_lp = p_ode[1]
    kp_pll = p_ode[2]
    ki_pll = p_ode[3]
    Ta = p_ode[4]
    kd = p_ode[5]
    kω = p_ode[6]
    kq = p_ode[7]
    ωf = p_ode[8]
    kpv = p_ode[9]
    kiv = p_ode[10]
    kffv = p_ode[11]
    rv = p_ode[12]
    lv = p_ode[13]
    kpc = p_ode[14]
    kic = p_ode[15]
    kffi = p_ode[16]
    ωad = p_ode[17]
    kad = p_ode[18]
    lf = p_ode[19]
    rf = p_ode[20]
    cf = p_ode[21]
    lg = p_ode[22]
    rg = p_ode[23]
    Vref = p_ode[24]    #Reference treated as parameter, but should NOT get updated during training
    ωref = p_ode[25]    #Reference treated as parameter, but should NOT get updated during training
    Pref = p_ode[26]    #Reference treated as parameter, but should NOT get updated during training
    Qref = p_ode[27]    #Reference treated as parameter, but should NOT get updated during training
    Xtrans = p_ode[28]
    Rtrans = p_ode[29]

    #STATE INDEX AND STATES
    i__vi_filter, vi_filter = 1, x[1]
    i__γd_ic, γd_ic = 2, x[2]
    i__vq_pll, vq_pll = 3, x[3]
    i__γq_ic, γq_ic = 4, x[4]
    i__ir_filter, ir_filter = 5, x[5]
    i__ξd_ic, ξd_ic = 6, x[6]
    i__ϕd_ic, ϕd_ic = 7, x[7]
    i__ε_pll, ε_pll = 8, x[8]
    i__ir_cnv, ir_cnv = 9, x[9]
    i__vr_filter, vr_filter = 10, x[10]
    i__ω_oc, ω_oc = 11, x[11]
    i__ξq_ic, ξq_ic = 12, x[12]
    i__vd_pll, vd_pll = 13, x[13]
    i__q_oc, q_oc = 14, x[14]
    i__ϕq_ic, ϕq_ic = 15, x[15]
    i__θ_pll, θ_pll = 16, x[16]
    i__θ_oc, θ_oc = 17, x[17]
    i__ii_cnv, ii_cnv = 18, x[18]
    i__ii_filter, ii_filter = 19, x[19]
    i__ir_nn, ir_nn = 20, x[20]
    i__ii_nn, ii_nn = 21, x[21]
    i__ir_offset, ir_offset = 22, x[22]
    i__ii_offset, ii_offset =  23, x[23]
    i__ir_out, ir_out =  24, x[24]
    i__ii_out, ii_out =  25, x[25]

    ω_base = 60.0*2*pi
    ω_sys = 1.0

    #PLL
    δω_pll = kp_pll * atan(vq_pll/vd_pll) + ki_pll* ε_pll
    ω_pll = δω_pll + ω_sys
    vd_filt_pll = sin(θ_pll + pi/2)*vr_filter - cos(θ_pll + pi/2)*vi_filter
    vq_filt_pll = cos(θ_pll + pi/2)*vr_filter + sin(θ_pll + pi/2)*vi_filter

    dx[i__vd_pll]=  ω_lp*(vd_filt_pll - vd_pll)                                 #docs:(1a)
    dx[i__vq_pll]=  ω_lp*(vq_filt_pll - vq_pll)                                 #docs:(1b)
    dx[i__ε_pll] =  atan(vq_pll/vd_pll)                                         #docs:(1c)
    dx[i__θ_pll] =  ω_base * δω_pll                                             #docs:(1d)

    #OUTER LOOP CONTROL
    pe = vr_filter*ir_filter + vi_filter*ii_filter                              #docs:(2d)
    qe = vi_filter*ir_filter - vr_filter*ii_filter                              #docs:(2e)
    v_ref_olc = Vref + kq * (Qref - q_oc)


    dx[i__ω_oc]= (Pref - pe - kd*(ω_oc - ω_pll) - kω*(ω_oc - ωref)) / Ta        #docs:(2a)
    dx[i__θ_oc]= ω_base*(ω_oc-ω_sys)                                            #docs:(2b)
    dx[i__q_oc]= ωf * (qe - q_oc)                                               #docs:(2c)

    #INNER LOOP CONTROL
    #reference transormations
    vd_filt_olc = sin(θ_oc + pi/2)*vr_filter - cos(θ_oc + pi/2)*vi_filter
    vq_filt_olc = cos(θ_oc + pi/2)*vr_filter + sin(θ_oc + pi/2)*vi_filter
    id_filt_olc = sin(θ_oc + pi/2)*ir_filter - cos(θ_oc + pi/2)*ii_filter
    iq_filt_olc = cos(θ_oc + pi/2)*ir_filter + sin(θ_oc + pi/2)*ii_filter
    id_cnv_olc =  sin(θ_oc + pi/2)*ir_cnv - cos(θ_oc + pi/2)*ii_cnv
    iq_cnv_olc =  cos(θ_oc + pi/2)*ir_cnv + sin(θ_oc + pi/2)*ii_cnv

    #Voltage control equations
    Vd_filter_ref = v_ref_olc - rv * id_filt_olc + ω_oc * lv * iq_filt_olc      #docs:(3g)
    Vq_filter_ref = -rv * iq_filt_olc - ω_oc * lv * id_filt_olc                 #docs:(3h)
    dx[i__ξd_ic]=  Vd_filter_ref -  vd_filt_olc                                 #docs:(3a)
    dx[i__ξq_ic]=  Vq_filter_ref -  vq_filt_olc                                 #docs:(3b)

    #current control equations
    Id_cnv_ref =(
      kpv*(Vd_filter_ref - vd_filt_olc) + kiv*ξd_ic -                           #docs:(3i)
      cf*ω_oc*vq_filt_olc + kffi*id_filt_olc
    )
    Iq_cnv_ref =(
      kpv*(Vq_filter_ref - vq_filt_olc) + kiv*ξq_ic +                           #docs:(3j)
      cf*ω_oc*vd_filt_olc + kffi*iq_filt_olc
    )
    dx[i__γd_ic]=  Id_cnv_ref - id_cnv_olc                                      #docs:(3c)
    dx[i__γq_ic]=  Iq_cnv_ref - iq_cnv_olc                                      #docs:(3d)

    #active damping equations
    Vd_cnv_ref =(
      kpc*(Id_cnv_ref - id_cnv_olc) + kic*γd_ic - lf*ω_oc*iq_cnv_olc +          #docs:(3k)
      kffv*vd_filt_olc - kad*(vd_filt_olc-ϕd_ic)
    )
    Vq_cnv_ref =(
      kpc*(Iq_cnv_ref - iq_cnv_olc) + kic*γq_ic + lf*ω_oc*id_cnv_olc +          #docs:(3l)
      kffv*vq_filt_olc - kad*(vq_filt_olc-ϕq_ic)
    )
    dx[i__ϕd_ic] =  ωad * (vd_filt_olc - ϕd_ic)                                 #docs:(3e)
    dx[i__ϕq_ic]=  ωad * (vq_filt_olc - ϕq_ic)                                  #docs:(3f)

    #LCL FILTER
    #reference transformations
    Vr_cnv =   sin(θ_oc + pi/2)*Vd_cnv_ref + cos(θ_oc + pi/2)*Vq_cnv_ref
    Vi_cnv =  -cos(θ_oc + pi/2)*Vd_cnv_ref + sin(θ_oc + pi/2)*Vq_cnv_ref

    Vr_pcc = Vm(t)*cos(Vθ(t)) + (ir_out*Rtrans - ii_out*Xtrans)
    Vi_pcc = Vm(t)*sin(Vθ(t)) + (ir_out*Xtrans + ii_out*Rtrans)

    dx[i__ir_cnv]= (ω_base/lf) *                                                #docs:(5a)
      (Vr_cnv - vr_filter - rf*ir_cnv + ω_sys*lf*ii_cnv)
    dx[i__ii_cnv]= (ω_base/lf) *                                                #docs:(5b)
      (Vi_cnv - vi_filter - rf*ii_cnv - ω_sys*lf*ir_cnv)
    dx[i__vr_filter]= (ω_base/cf) *                                             #docs:(5c)
      (ir_cnv - ir_filter + ω_sys*cf*vi_filter)
    dx[i__vi_filter]= (ω_base/cf) *                                             #docs:(5d)
      (ii_cnv - ii_filter - ω_sys*cf*vr_filter)
    dx[i__ir_filter]= (ω_base/lg) *                                             #docs:(5e)
      (vr_filter - Vr_pcc - rg*ir_filter + ω_sys*lg*ii_filter)
    dx[i__ii_filter]= +(ω_base/lg) *                                            #docs:(5f)
      (vi_filter - Vi_pcc - rg*ii_filter - ω_sys*lg*ir_filter)

     #NN CURRENT SOURCE (NN input includes all gfm states)
     dx[i__ir_nn] = ir_offset -  nn_states(x[1:19],p_nn)[1]
     dx[i__ii_nn] = ii_offset -  nn_states(x[1:19],p_nn)[2]

     #Current offset for the NN Current Source
     dx[i__ir_offset] = 0.0
     dx[i__ii_offset] = 0.0

     #ALEGRAIC STATE - OUTPUT CURRENT
     dx[i__ir_out] =  ir_out - ir_nn - ir_filter
     dx[i__ii_out] =  ii_out - ii_nn - ii_filter
end
