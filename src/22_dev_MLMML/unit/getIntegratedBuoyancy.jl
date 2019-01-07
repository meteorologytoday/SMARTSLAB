using Test
include("../MLMML.jl")
using .MLMML

zs = collect(Float64, range(0, stop=-1000, length=101))
bs = zeros(Float64, length(zs)-1)
h = 133.0
b_ML = 10.0
Δb   = 7.0
L = 500.0
for i = 1:length(bs)
    if zs[i] > -h
        bs[i] = b_ML
    else
        z = (zs[i] + zs[i+1]) / 2.0
        bs[i] = (b_ML - Δb) * exp((z+h)/L)
    end
end

test_zs = collect(Float64, range(0, stop=-1000, length=2001))
test_ans = copy(test_zs)

test_res = copy(test_zs)
for i = 1:length(test_ans)
    if test_zs[i] > -h
        test_ans[i] = b_ML * (-test_zs[i])
    else
        test_ans[i] = b_ML * h + (b_ML - Δb) * L * (1.0 - exp((test_zs[i]+h)/L))
    end

    test_res[i] = MLMML.getIntegratedBuoyancy(zs=zs, bs=bs, b_ML=b_ML, h=h, target_z=test_zs[i])
end


using PyPlot
fig, ax = plt[:subplots](1, 2, figsize=(8,6), sharey=true)

ax[1][:plot](bs, (zs[1:end-1] + zs[2:end])/2, "k-", label="bs")
ax[1][:legend]()

ax[2][:plot](test_ans, test_zs, "r--", label="test_answer")
ax[2][:plot](test_res, test_zs, "b-.", label="test_result")
ax[2][:legend]()
plt[:show]()


# Test begin
println("Test Begin")
for i = 1:length(test_ans)
    I = MLMML.getIntegratedBuoyancy(zs=zs, bs=bs, b_ML=b_ML, h=h, target_z=test_zs[i])
    @test I ≈ test_ans[i] rtol=.05
end
println("Test Done")

