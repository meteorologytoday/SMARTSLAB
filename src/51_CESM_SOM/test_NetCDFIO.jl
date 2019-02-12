include("julia_lib/NetCDFIO.jl")

using .NetCDFIO

# 10 deg
domain_file = "/home/tienyiah/projects/cesm2_test/inputdata/share/domains/domain.ocn.fv10x15_gx3v7.180321.nc"

# 2 deg
#domain_file = "~/projects/cesm2_test/inputdata/share/domains/domain.ocn.1.9x2.5_gx1v6_090403.nc"

output_file = "test.nc"

mi = NetCDFIO.MapInfo{Float64}(domain_file)
NetCDFIO.createNCFile(mi, output_file)


sst = rand(size(mi.xc)..., 10) * 50
mld = rand(size(mi.xc)..., 10) * 50

#NetCDFIO.write2NCFile(mi, output_file, "sst", sst)

for i = 1:300
    sst = rand(size(mi.xc)...) * 3 .+ i
    mld = rand(size(mi.xc)...) * 3 .+ i^2
    NetCDFIO.write2NCFile(mi, output_file, "sst", sst; time=i)
    NetCDFIO.write2NCFile(mi, output_file, "mld", mld; time=i)
end





# fake some sst data
