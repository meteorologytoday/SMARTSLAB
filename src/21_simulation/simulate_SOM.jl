include("config.jl")
include("../lib/simulate_cores/SOM.jl")

using NCDatasets
using .KTSimulation

hQ_nc_filename = "HMC_SST_Td-fixed_NCAR_5deg_init-30m_c4_s1000_w200.nc" 
hQ_nc_filename = normpath(joinpath(data_path, hQ_nc_filename)) 

sim_nc_filename = format("21_KTsimulate_{}.nc", splitext(basename(hQ_nc_filename))[1])
sim_nc_filename = joinpath(data_path, sim_nc_filename)


sim_len = 1200
θd = 273.15 * ρ * c_p


SST = zeros(dtype, length(lon), length(lat), sim_time)
SST .= NaN

ds = Dataset(hQ_nc_filename,"r")
h = nomissing(ds["h_mean"][:], NaN)
Q = nomissing(ds["Q_mean"][:], NaN)
close(ds)

θ = readModelVar("tos") * ρ * c_p
F = readModelVar("hfds")



println("Size of h / Q : ", size(h))

for i = 1:length(lon), j = 1:length(lat)
    if isnan(h[i, j, 1])
        continue
    end

    data = KTSimulation.run(;
        θ_init  = θ[i, j, init_time],
        θd      = θd,
        Δt      = Δt,
        F       = F[i, j, init_time:init_time + total_sim_time - 1],
        Q       = Q[i, j, :],
        h       = h[i, j, :],
        period  = 12,
        ret_len = total_sim_time,
    )

    SST[i, j, :] = data["θ"][warmup_time + 1 : warmup_time + sim_time]

end

SST /= (ρ * c_p)



nan2missing!(SST)

ds = Dataset(sim_nc_filename,"c")
defDim(ds,"time", sim_len)
defDim(ds,"lat", length(lat))
defDim(ds,"lon", length(lon))

defVar(ds, "time", Float64, ("time",))[:] = collect(1:sim_len)
defVar(ds, "lat", Float64, ("lat",))[:] = lat
defVar(ds, "lon", Float64, ("lon",))[:] = lon

for o in (
    [
        "SST", SST, ("lon", "lat", "time"), Dict(
        "long_name"=>"SST predicted",
        "units" => "K",
        )
    ], 
)
    varname, vardata, vardims, varatts = o
    println("Writing ", varname, " with size: ", size(vardata), " ; dim: ", vardims)

    ncvar = defVar(ds, varname, eltype(vardata), vardims)
    ncvar.attrib["_FillValue"] = missing_value
    for key in keys(varatts)
        ncvar.attrib[key] = varatts[key]
    end

    ncvar[:] = vardata
    println("done.")
end

close(ds)
println("Output file: ", sim_nc_filename)
