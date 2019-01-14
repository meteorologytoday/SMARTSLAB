function doDiffusion!(;
    hs   :: Array{Float64, 1}, 
    zs   :: Array{Float64, 1}, 
    Δzs  :: Array{Float64, 1}, 
    bs   :: Array{Float64, 1}, 
    b_ML :: Float64,
    h    :: Float64,
    FLDO :: Integer,
)

    n = length(bs) - (FLDO - 1)
    A = spzeros(Float64, n, n)
    

    



end
