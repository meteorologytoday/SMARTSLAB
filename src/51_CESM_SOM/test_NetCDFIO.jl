include("julia_lib/NetCDFIO.jl")

using .NetCDFIO

# 10 deg
domain_file = "/home/tienyiah/projects/cesm2_test/inputdata/share/domains/domain.ocn.fv10x15_gx3v7.180321.nc"

# 2 deg
#domain_file = "~/projects/cesm2_test/inputdata/share/domains/domain.ocn.1.9x2.5_gx1v6_090403.nc"

output_file = "test.nc"

mi = NetCDFIO.MapInfo{Float64}(domain_file)
NetCDFIO.createNCFile(mi, output_file)


sst = rand(size(mi.xc)..., 100) * 50
NetCDFIO.write2NCFile(mi, output_file, "sst", sst)

for i = 1:20
    sst = rand(size(mi.xc)...) * 50 .+ 100
    NetCDFIO.write2NCFile(mi, output_file, "sst", sst)
end


#sst = rand(size(mi.xc)..., 1) * 0.0
#NetCDFIO.write2NCFile(mi, output_file, "sst", sst)



# fake some sst data
