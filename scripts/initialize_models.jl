using NLsolve

f = get_init_gfm(p, vr_filter, vi_filter, ir_filter, ii_filter)

x0_init =   [vi_filter, # Same as above
            γd_ic,
            vq_pll,
            γq_ic,
            ir_filter,  # Same as above
            ξd_ic,
            ϕd_ic,
            ε_pll,
            ir_cnv,
            vr_filter,# Same as above
            ω_oc,
            ξq_ic,
            vd_pll,
            q_oc,
            ϕq_ic,
            θ_pll,
            θ_oc,
            ii_cnv,
            ii_filter# Same as above
        ]

nlsolve(f, j!, x0_init)
