include("single_point.jl")

function extend(a::AbstractArray, len::Int)
    if length(a) < len
        a = repeat(a, outer=ceil(Int, len / length(a)))
    end

    return a[1:len]
end


ret_len = 24
P = 1.0
K = 4.0 / P

t = collect(linspace(0.0, P, ret_len+1))[1:end-1]
println(t)

dt = t[2] - t[1]

S = t * 0.0
B = t * 0.0 + 1.0
G = t * 0.0 + 1.0
D = t * 0.0

for i = 1:ret_len
    _t = t[i]
    if _t < P / 4.0
        S[i] = K * _t
    elseif _t < P * 3.0 / 4.0
        S[i] = K * (P / 2.0 - _t)
    else
        S[i] = K * (_t - P)
    end
end


dSdt = extend(S, 3*ret_len)[ret_len:2*ret_len+1]
dSdt = (dSdt[3:end] - dSdt[1:end-2]) / dt



h_init  = 30.0
Ts_init = 273.15 + 10
Td      = 273.15
β       = Inf

h, Ts = KrausTurnerMLD_simplified(
    h_init,
    Ts_init,
    Td,
    dt,
    β,
    S,
    dSdt,
    B,
    G,
    D,
    ret_len 
)

println(h)
println(Ts)

using PyPlot

fig, ax = plt[:subplots](3, 1, sharex=true, figsize=(12, 8))

ax1_twin = ax[1][:twinx]()
ax2_twin = ax[2][:twinx]()

ax[1][:plot](t, S)
ax1_twin[:plot](t, dSdt)

ax[2][:plot](t, Ts)
ax2_twin[:plot](t, h)



plt[:show]()
