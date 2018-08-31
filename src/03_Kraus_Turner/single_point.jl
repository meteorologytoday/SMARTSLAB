


function f(
    h       :: TYPE,
    Ts      :: TYPE,
    Td      :: TYPE,
    β       :: TYPE,
    S       :: TYPE,
    dSdt    :: TYPE,
    B       :: TYPE,
    G       :: TYPE,
    D       :: TYPE,
) where TYPE <: AbstractFloat

    Γ = G - D + S / β
    Q = S + B

    diag_h = nothing
    
    # test if dhdt > 0
    dhdt = 1.0 / Q * (2.0 / β - h) * dSdt

    if dhdt <= 0 # h can be diagnosed
        dTsdt = 0.5 * Q^2.0 / Γ
        diag_h = 2.0 * Γ / Q
    else
        dTsdt = 2.0 / h^2.0 * (Q * h - Γ)
        dhdt  = 1.0 / (Ts - Td) / h * (2.0 * Γ - Q * h)
    end


    return diag_h, dhdt, dTsdt

end






function KrausTurnerMLD_simplified(
    h_init  :: TYPE,
    Ts_init :: TYPE,
    Td      :: TYPE,
    dt      :: TYPE,
    β       :: TYPE,
    S       :: Array{TYPE, 1},
    dSdt    :: Array{TYPE, 1},
    B       :: Array{TYPE, 1},
    G       :: Array{TYPE, 1},
    D       :: Array{TYPE, 1},
    ret_len :: Integer
) where TYPE <: AbstractFloat

    S    = extend(S,    ret_len)
    B    = extend(B,    ret_len)
    G    = extend(G,    ret_len)
    D    = extend(D,    ret_len)
    dSdt = extend(dSdt, ret_len)


    h  = zeros(TYPE, ret_len)
    Ts = h * 0.0

    h[1]  = h_init
    Ts[1] = Ts_init

    for i = 1 : ret_len-1


        h_diag, k1_dhdt, k1_dTsdt = f(
            h[i],
            Ts[i],
            Td,
            β,
            S[i],
            dSdt[i],
            B[i],
            G[i],
            D[i]
        )

        h[i+1]  = h[i]  + k1_dhdt  * dt
        Ts[i+1] = Ts[i] + k1_dTsdt * dt


    end


    return h, Ts
end




