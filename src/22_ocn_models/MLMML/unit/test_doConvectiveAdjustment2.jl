using Test
include("../MLMML.jl")
using .MLMML

zs = collect(Float64, range(0, stop=-1000, length=2001))
h = MLMML.h_min
b_ML = 10.0
Δb   = -1e-2
L = 500.0

oc = MLMML. makeSimpleOceanColumn(zs=zs)
oc.bs[oc.FLDO:end] .-= Δb

oc_1 = MLMML.copy(oc)
MLMML.doConvectiveAdjustment!(oc)
oc_2 = MLMML.copy(oc)


test_zs = collect(Float64, range(0, stop=-1000, length=2001))
test_intb_1 = test_zs * 0.0
test_intb_2 = test_zs * 0.0
for i = 1:length(test_zs)
    test_intb_1[i] = MLMML.getIntegratedBuoyancy(oc_1; target_z=test_zs[i])
    test_intb_2[i] = MLMML.getIntegratedBuoyancy(oc_2; target_z=test_zs[i])
end

using PyPlot
fig, ax = plt[:subplots](1, 3, figsize=(8,6), sharey=true)

ax[1][:plot](oc_1.bs, (zs[1:end-1] + zs[2:end])/2, "k-", label="bs_1")
ax[1][:plot](oc_2.bs, (zs[1:end-1] + zs[2:end])/2, "r--", label="bs_2")
ax[1][:legend]()

ax[2][:plot](test_intb_1, test_zs, "k-", label="intb_1")
ax[2][:plot](test_intb_2, test_zs, "r--", label="intb_2")
ax[2][:legend]()

ax[3][:plot](test_intb_2 - test_intb_1, test_zs, "b", label="intb_2 - intb_1")
ax[3][:legend]()

plt[:show]()

