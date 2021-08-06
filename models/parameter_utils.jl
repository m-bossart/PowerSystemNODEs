
"""
    randomize_parameters!(inv::DynamicInverter, param_range::Tuple{Float64,Float64})

Change each parameter by scaling by randomly generated factor in param_range
"""
function randomize_parameters!(inv::DynamicInverter, param_range::Tuple{Float64,Float64})
    #pll
    set_ω_lp!(inv.freq_estimator, get_ω_lp(inv.freq_estimator)*rand(Uniform(param_range[1],param_range[2])))
    set_kp_pll!(inv.freq_estimator, get_kp_pll(inv.freq_estimator)*rand(Uniform(param_range[1],param_range[2])))
    set_ki_pll!(inv.freq_estimator, get_ki_pll(inv.freq_estimator)*rand(Uniform(param_range[1],param_range[2])))
    #outer control
    PSY.set_Ta!(inv.outer_control.active_power, PSY.get_Ta(inv.outer_control.active_power)*rand(Uniform(param_range[1],param_range[2])))
    PSY.set_kd!(inv.outer_control.active_power, PSY.get_kd(inv.outer_control.active_power)*rand(Uniform(param_range[1],param_range[2])))
    PSY.set_kω!(inv.outer_control.active_power, PSY.get_kω(inv.outer_control.active_power)*rand(Uniform(param_range[1],param_range[2])))
    PSY.set_kq!(inv.outer_control.reactive_power, PSY.get_kq(inv.outer_control.reactive_power)*rand(Uniform(param_range[1],param_range[2])))
    PSY.set_ωf!(inv.outer_control.reactive_power, PSY.get_ωf(inv.outer_control.reactive_power)*rand(Uniform(param_range[1],param_range[2])))
    #inner control
    set_kpv!(inv.inner_control, get_kpv(inv.inner_control)*rand(Uniform(param_range[1],param_range[2])))
    set_kiv!(inv.inner_control, get_kiv(inv.inner_control)*rand(Uniform(param_range[1],param_range[2])))
    set_kffv!(inv.inner_control, get_kffv(inv.inner_control)*rand(Uniform(param_range[1],param_range[2])))
    set_rv!(inv.inner_control, get_rv(inv.inner_control)*rand(Uniform(param_range[1],param_range[2])))
    set_lv!(inv.inner_control, get_lv(inv.inner_control)*rand(Uniform(param_range[1],param_range[2])))
    set_kpc!(inv.inner_control, get_kpc(inv.inner_control)*rand(Uniform(param_range[1],param_range[2])))
    set_kic!(inv.inner_control, get_kic(inv.inner_control)*rand(Uniform(param_range[1],param_range[2])))
    set_kffi!(inv.inner_control, get_kffi(inv.inner_control)*rand(Uniform(param_range[1],param_range[2])))
    set_ωad!(inv.inner_control, get_ωad(inv.inner_control)*rand(Uniform(param_range[1],param_range[2])))
    set_kad!(inv.inner_control, get_kad(inv.inner_control)*rand(Uniform(param_range[1],param_range[2])))
    #lcl
    set_lf!(inv.filter, get_lf(inv.filter)*rand(Uniform(param_range[1],param_range[2])))
    set_rf!(inv.filter, get_rf(inv.filter)*rand(Uniform(param_range[1],param_range[2])))
    set_cf!(inv.filter, get_cf(inv.filter)*rand(Uniform(param_range[1],param_range[2])))
    set_lg!(inv.filter, get_lg(inv.filter)*rand(Uniform(param_range[1],param_range[2])))
    set_rg!(inv.filter, get_rg(inv.filter)*rand(Uniform(param_range[1],param_range[2])))
end

function set_parameters!(inv::DynamicInverter, p::Vector{Float64})
    #pll
    set_ω_lp!(inv.freq_estimator, p[1])
    set_kp_pll!(inv.freq_estimator, p[2])
    set_ki_pll!(inv.freq_estimator ,p[3])
    #outer control
    PSY.set_Ta!(inv.outer_control.active_power, p[4])
    PSY.set_kd!(inv.outer_control.active_power, p[5])
    PSY.set_kω!(inv.outer_control.active_power, p[6])
    PSY.set_kq!(inv.outer_control.reactive_power, p[7])
    PSY.set_ωf!(inv.outer_control.reactive_power,p[8])
    #inner control
    set_kpv!(inv.inner_control, p[9])
    set_kiv!(inv.inner_control, p[10])
    set_kffv!(inv.inner_control, p[11])
    set_rv!(inv.inner_control, p[12])
    set_lv!(inv.inner_control, p[13])
    set_kpc!(inv.inner_control, p[14])
    set_kic!(inv.inner_control, p[15])
    set_kffi!(inv.inner_control, p[16])
    set_ωad!(inv.inner_control, p[17])
    set_kad!(inv.inner_control, p[18])
    #lcl
    set_lf!(inv.filter, p[19])
    set_rf!(inv.filter, p[20])
    set_cf!(inv.filter, p[21])
    set_lg!(inv.filter, p[22])
    set_rg!(inv.filter, p[23])
end

function get_parameters(inv::DynamicInverter)
    #pll
    p = Vector{Float64}(undef,23)
    p[1] = get_ω_lp(inv.freq_estimator)
    p[2] = get_kp_pll(inv.freq_estimator)
    p[3] = get_ki_pll(inv.freq_estimator)
    #outer control
    p[4] = PSY.get_Ta(inv.outer_control.active_power)
    p[5] = PSY.get_kd(inv.outer_control.active_power)
    p[6] = PSY.get_kω(inv.outer_control.active_power)
    p[7] = PSY.get_kq(inv.outer_control.reactive_power)
    p[8] = PSY.get_ωf(inv.outer_control.reactive_power)
    #inner control
    p[9] = get_kpv(inv.inner_control)
    p[10] = get_kiv(inv.inner_control)
    p[11] = get_kffv(inv.inner_control)
    p[12] = get_rv(inv.inner_control)
    p[13] = get_lv(inv.inner_control)
    p[14] = get_kpc(inv.inner_control)
    p[15] = get_kic(inv.inner_control)
    p[16] = get_kffi(inv.inner_control)
    p[17] = get_ωad(inv.inner_control)
    p[18] = get_kad(inv.inner_control)
    #lcl
    p[19] = get_lf(inv.filter)
    p[20] = get_rf(inv.filter)
    p[21] = get_cf(inv.filter)
    p[22] = get_lg(inv.filter)
    p[23] = get_rg(inv.filter)
    return p
end
