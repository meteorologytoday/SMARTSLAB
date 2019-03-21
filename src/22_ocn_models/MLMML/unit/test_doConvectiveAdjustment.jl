using Test
include("../MLMML.jl")
using .MLMML

zs = collect(Float64, range(0, stop=-1000, length=2001))
ocs = []
adj_ocs = []

# OC1
oc1 = MLMML.makeSimpleOceanColumn(zs=zs, Δb=-5e-3)
push!(ocs, oc1)

# OC2
oc2 = MLMML.makeSimpleOceanColumn(zs=zs, b_slope = -30.0 / 5000.0 * MLMML.g * MLMML.α)
push!(ocs, oc2)

# OC3
k  = 2π / 500.0 * 3
amp = 1e-3
oc3 = MLMML.copy(oc1)

mid_zs = (oc3.zs[1:end-1] + oc3.zs[2:end]) / 2.0
oc3.bs[oc3.FLDO:end] += amp * sin.( (mid_zs[oc3.FLDO:end] .- mid_zs[oc3.FLDO]) * k )
push!(ocs, oc3)

# OC4
noise_amp = amp / 2.0
oc4 = MLMML.copy(oc3)
oc4.bs[oc4.FLDO:end] += noise_amp * randn(length(oc4.bs) - oc4.FLDO + 1)
push!(ocs, oc4)


for i = 1:length(ocs)
    oc_copy = MLMML.copy(ocs[i])
    MLMML.OC_doConvectiveAdjustment!(oc_copy)
    push!(adj_ocs, oc_copy)
end

function genLineCoord(oc::MLMML.OceanColumn)

    x = []
    z = []

    push!(x, oc.b_ML)
    push!(z, 0.0, - oc.h)

    if oc.FLDO != -1 && oc.FLDO < length(oc.bs)
        push!(x, oc.bs[oc.FLDO])
        push!(z, -oc.h, oc.zs[oc.FLDO+1])
        for i=oc.FLDO+1:length(oc.bs)
            push!(x, oc.bs[i])
            push!(z, oc.zs[i], oc.zs[i+1])
        end
    end

    x = repeat(x, inner=(2,))

    return x, z
end


using PyPlot
using Formatting

for i = 1:length(ocs)

    oc = ocs[i]
    adj_oc = adj_ocs[i]

    fig, ax = plt[:subplots](1, 1, figsize=(8,6), sharey=true)

    ax[:plot](genLineCoord(oc)..., "k-", label="Original")
    ax[:plot](genLineCoord(adj_oc)..., "r--", label="Adjusted")

    ax[:set_title](
        format("Integrated buoyancy: before={:.2e}, change={:.2e}",
            MLMML.getIntegratedBuoyancy(oc),
            MLMML.getIntegratedBuoyancy(adj_oc) - MLMML.getIntegratedBuoyancy(oc),
        )
    )
    
    
    
    ax[:legend]()

    plt[:show]()

end
