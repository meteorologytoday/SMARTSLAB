include("single_point.jl")


ret_len = 12 * 10






h, Ts = KrausTurnerMLD_simplified(
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
)
