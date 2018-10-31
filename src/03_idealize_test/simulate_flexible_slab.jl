include("core.jl")
using Printf
function prtArr(A)
    for i = 1:size(A)[1]
        for j = 1:size(A)[2]

            @printf("%.1f ", A[i,j]) 

        end
        @printf("\n")
    end
end


ρ    = 1027.0  # kg / m^3
c_p  = 3985.0  # J / kg / K

dtype = typeof(ρ)


years = 10
months_per_year = 12
pts_per_month = 10
pts_per_year = months_per_year * pts_per_month
secs_per_year = 86400.0 * 365.0

ret_len = months_per_year * pts_per_month * years

t = collect(range(0.0, stop=years * secs_per_year, length=ret_len+1))[1:end-1]
Δt = t[2] - t[1]

t_year = t / secs_per_year

S = zeros(dtype, pts_per_year)
h = copy(S)
B = copy(S)

S0 = 100.0
for i = 1 : length(S)
    S[i] = S0 * sin(2π * t_year[i])
    h[i] = 30.0 - 10.0 * sin(2π * t_year[i])
end

h .= 30.0

S    = extend(S,    ret_len)
B    = extend(B,    ret_len)
h    = extend(h,    ret_len)
 
Ts_init = 300.00 
Td      = 0.0
β       = Inf


θs_init = Ts_init * ρ * c_p
θd = Td * ρ * c_p

θs = simulation(;
    θs_init = θs_init,
    θd      = θd,
    Δt      = Δt,
    β       = β,
    S       = S,
    B       = B,
    h       = h,
    ret_len = ret_len
)

Ts = θs / ρ / c_p

function avg_every_pts(
    x :: Array{T},
    pts :: Int
) where T <: AbstractFloat
    N = length(x)

    if mod(N, pts) != 0
        throw(Exception("Length of array must be multiple of pts input."))
    end

    xx = zeros(eltype(x), convert(Int, N / pts))
    
    for i = 1 : length(xx)
        xx[i] = sum(x[(i-1)*pts+1 : i*pts]) / pts
    end

    return xx
end

t = avg_every_pts(t, pts_per_month)
t_year = avg_every_pts(t_year, pts_per_month)
Ts = avg_every_pts(Ts, pts_per_month)
θs = avg_every_pts(θs, pts_per_month)
h = avg_every_pts(h, pts_per_month)
S = avg_every_pts(S, pts_per_month)
B = avg_every_pts(B, pts_per_month)
Δt = t[2] - t[1]

A = S0 * secs_per_year / (2π * 30.0)
c0 = A + θs_init
Ts_analytic = (- A * cos.(2π * t_year) .+ c0) / ρ / c_p


println("True h: ", h)
println("Δt:  ", Δt)



include("../lib/LR_shallow_water.jl")

LR_h, LR_Q = LR_shallow_water(
    years = years,
    pts_per_year = months_per_year,
    F = S + B,
    θ = θs .- θd,
    Δt = Δt 
)


LR_h = extend(LR_h, length(t_year))



















