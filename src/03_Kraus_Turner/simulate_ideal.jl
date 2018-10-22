using Printf

include("single_point.jl")

function extend(a::AbstractArray, len::Int)
    if length(a) < len
        a = repeat(a, outer=ceil(Int, len / length(a)))
    end

    return a[1:len]
end

ρ    = 1027.0  # kg / m^3
c_p  = 3985.0  # J / kg / K

MAXQ = 344.0
ret_len = 100
P = 1.0
K = 4.0 / P * MAXQ

t = collect(range(0.0, stop=P, length=ret_len+1))[1:end-1]

dt = t[2] - t[1]

S = t * 0.0
B = t * 0.0 .+ 1.0
G = t * 0.0 .+ 1.0
D = t * 0.0

for i = 1:ret_len
    _t = t[i]
    if _t < P / 4.0
        S[i] = K * _t
    elseif _t < P * 3.0 / 4.0
        S[i] = K * (P / 2.0 .- _t)
    else
        S[i] = K * (_t .- P)
    end
end

S = circshift(S, floor(ret_len/3))

dSdt = (circshift(S, -1) - circshift(S, 1) ) / 2.0 / dt
dBdt = (circshift(B, -1) - circshift(B, 1) ) / 2.0 / dt
dGdt = (circshift(G, -1) - circshift(G, 1) ) / 2.0 / dt
dDdt = (circshift(D, -1) - circshift(D, 1) ) / 2.0 / dt


@printf("dt: %f\n", dt)


h_init  = 30.0
Ts_init = 273.15 + 30
Td      = 273.15
β       = 1000000.0

θs_init = Ts_init * c_p * ρ
θd      = Td * c_p * ρ


h, θs = KrausTurnerMLD_simplified(
    h_init = h_init,
    θs_init = θs_init,
    θd = θd,
    dt = dt,
    β = β,
    S = S,
    B = B,
    G = G,
    D = D,
    dSdt = dSdt,
    dBdt = dBdt,
    dGdt = dGdt,
    dDdt = dDdt,
    ret_len = ret_len
)

Ts = θs / c_p / ρ

println(h)
println(Ts)

using PyPlot

fig, ax = plt[:subplots](3, 1, sharex=true, figsize=(12, 8))

ax1_twin = ax[1][:twinx]()
ax2_twin = ax[2][:twinx]()

ax[1][:set_title]("S and dSdt")
ax[1][:plot](t, S, color="k", label="S")
ax1_twin[:plot](t, dSdt, color="k", dashes=(5,3), label="dSdt")
ax[1][:legend]()
ax1_twin[:legend]()


ax[2][:plot](t, Ts, color="k", label="Ts")
ax[2][:legend]()

ax[3][:plot](t, -h, color="k", label="h")
ax[3][:legend]()






plt[:show]()
