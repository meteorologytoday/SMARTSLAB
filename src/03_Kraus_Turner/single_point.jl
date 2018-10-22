


function f(;
    h       :: TYPE,
    θs      :: TYPE,
    θd      :: TYPE,
    β       :: TYPE,
    S       :: TYPE,
    B       :: TYPE,
    G       :: TYPE,
    D       :: TYPE,
    dSdt    :: TYPE,
    dBdt    :: TYPE,
    dGdt    :: TYPE,
    dDdt    :: TYPE,
) where TYPE <: AbstractFloat

    Γ    = G - D + S / β
    dΓdt = dGdt - dDdt + dSdt / β
    Q    = S + B
    dQdt = dSdt + dBdt

    diag_h = nothing
    
    # test if dhdt > 0
    #dhdt_assume = 1.0 / Q * (2.0 / β - h) * dSdt
    
    dhdt_assume = 2.0 * ((dΓdt * Q - Γ * dQdt) / Q^2.0)

    dhdt = dhdt_assume
    if dhdt_assume <= 0 # h can be diagnosed
        println("Shallowing")
        dθsdt = 0.5 * Q^2.0 / Γ
        diag_h = 2.0 * Γ / Q
    else
        println("Deepening")
        dθsdt = 2.0 / h^2.0 * (Q * h - Γ)
        dhdt  = 1.0 / (θs - θd) / h * (2.0 * Γ - Q * h)
    end

    @printf("dhdt_assume: %f, dhdt: %f, dθsdt: %f, Γ: %f, S: %f, Q: %f, dSdt: %f\n", dhdt_assume, dhdt, dθsdt, Γ, S, Q, dSdt)



    return diag_h, dhdt, dθsdt

end






function KrausTurnerMLD_simplified(;
    h_init  :: TYPE,
    θs_init :: TYPE,
    θd      :: TYPE,
    dt      :: TYPE,
    β       :: TYPE,
    S       :: Array{TYPE, 1},
    B       :: Array{TYPE, 1},
    G       :: Array{TYPE, 1},
    D       :: Array{TYPE, 1},
    dSdt    :: Array{TYPE, 1},
    dBdt    :: Array{TYPE, 1},
    dGdt    :: Array{TYPE, 1},
    dDdt    :: Array{TYPE, 1},
    ret_len :: Integer
) where TYPE <: AbstractFloat

    S    = extend(S,    ret_len)
    B    = extend(B,    ret_len)
    G    = extend(G,    ret_len)
    D    = extend(D,    ret_len)
    dSdt = extend(dSdt, ret_len)
    dBdt = extend(dBdt, ret_len)
    dGdt = extend(dGdt, ret_len)
    dDdt = extend(dDdt, ret_len)


    h  = zeros(TYPE, ret_len)
    θs = copy(h)

    h[1]  = h_init
    θs[1] = θs_init

    for i = 1 : ret_len-1
        
        h_diag, k1_dhdt, k1_dθsdt = f(
            h=h[i],
            θs=θs[i],
            θd=θd,
            β=β,
            S=S[i],
            B=B[i],
            G=G[i],
            D=D[i],
            dSdt=dSdt[i],
            dBdt=dBdt[i],
            dGdt=dGdt[i],
            dDdt=dDdt[i],
        )

        h[i+1]  = h[i]  + k1_dhdt  * dt
        θs[i+1] = θs[i] + k1_dθsdt * dt

    end


    return h, θs
end




