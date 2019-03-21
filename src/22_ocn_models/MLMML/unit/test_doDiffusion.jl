using Test
include("../MLMML.jl")
using .MLMML

zs = collect(Float64, range(0, stop=-100, length=53))

ts = collect(Float64, range(0, step=86400.0, length=2))
Δt = ts[2] - ts[1]

#oc = MLMML.makeBlankOceanColumn(zs=zs)

#oc.K = 1e-5
#oc.bs[:] .= 0.0
#=
MLMML.setMixedLayer!(
    oc;
    b_ML = 1.0,
    h = 30.0
)
=#

Δb_init = 0.5 * 10.0 * MLMML.α * MLMML.g  
b_ML_init = 1.0
h_init = 50.0
K = 1e-5

oc = MLMML.makeSimpleOceanColumn(
    zs      = zs,
    b_slope = b_slope,
    b_ML    = b_ML_init,
    h       = h_init,
    Δb      = Δb_init,
    K       = K
)

oc.bs[:] .= 0.0
oc.bs[20:30] .= 1.0

oc.bs += oc.zs[1:end-1] * 10.0 * MLMML.g * MLMML.α 
MLMML.setMixedLayer!(
    oc;
    b_ML = 0.0,
    h = MLMML.h_min,
    convective_adjustment=false
)



ibuo     = zeros(length(ts))
bs_rec   = zeros(length(oc.bs), length(ts))
b_ML_rec = zeros(length(ts))
h_rec    = zeros(length(ts))
FLDO_rec = zeros(Integer, length(ts))

for i=1:length(ts)

    if i!=1
        MLMML.doDiffusion_EulerBackward!(oc; Δt=Δt)
    end
    ibuo[i] = MLMML.getIntegratedBuoyancy(oc) 
    bs_rec[:, i] = oc.bs
    b_ML_rec[i] = oc.b_ML
    h_rec[i] = oc.h
    FLDO_rec[i] = oc.FLDO
end



using PyPlot
fig, ax = plt[:subplots](1, 2, figsize=(12,6), squeeze=false)

function plotProfile(;
    h       :: Float64,
    b_ML    :: Float64,
    bs      :: Array{Float64,1},
    zs      :: Array{Float64,1},
    FLDO    :: Integer,
    linefmt :: String="k-"
)
    ax[1][:plot]([b_ML, b_ML], [0, -h], linefmt)
    ax[1][:plot]([bs[FLDO], bs[FLDO]], [-h, zs[FLDO+1]], linefmt)

    if FLDO < length(bs)
        for i=FLDO+1:length(bs)
            ax[1][:plot]([bs[i], bs[i]], [zs[i], zs[i+1]], linefmt)
        end
    end
end

plotProfile(
    h=h_rec[1],
    b_ML=b_ML_rec[1],
    FLDO=FLDO_rec[1],
    bs=bs_rec[:, 1],
    zs=zs,
    linefmt="k--"
)

for i=2:length(ts)
    plotProfile(
        h=h_rec[i],
        b_ML=b_ML_rec[i],
        FLDO=FLDO_rec[i],
        bs=bs_rec[:, i],
        zs=zs,
        linefmt="-"
    )
end

ax[2][:set_title]("Integrated buoyancy with time (should be a constant)")
ax[2][:plot](ts, ibuo, label="integrated buoyancy")
ax[2][:legend]()

plt[:show]()

