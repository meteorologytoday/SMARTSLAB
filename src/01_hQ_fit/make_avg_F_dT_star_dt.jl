include("config.jl")
include("NetCDFHelper.jl")

using NetCDF

dt2 = 2.0 * dt

T_star     = SST[:,:,13:end-12] * ρ * c_p
dT_star_dt = (SST[:,:,14:end-11] - SST[:,:,12:end-13]) * ρ * c_p / dt2
TOT_F      = TOT_F[:,:,13:end-12]

dT_star_dt_mean = zeros(eltype(TOT_F), size(TOT_F)[1:2]..., 12)
dT_star_dt_std  = copy(dT_star_dt_mean)

TOT_F_mean  = copy(dT_star_dt_mean)
TOT_F_std   = copy(dT_star_dt_mean)


T_star_mean  = copy(dT_star_dt_mean)

for m = 1:12
    @printf("Doing Month %02d\n", m)

    tup                      = (dT_star_dt[:, :, m:12:end], 3)
    dT_star_dt_mean[:, :, m] = mean(tup...)
    dT_star_dt_std[:, :, m]  = std(tup...)


    tup                      = (TOT_F[:, :, m:12:end], 3)
    TOT_F_mean[:, :, m]      = mean(tup...)
    TOT_F_std[:,:,m]         = std(tup...)
    
    T_star_mean[:, :, m]      = mean(T_star[:,:,m:12:end], 3)
end

mask = isnan.(dT_star_dt_mean)

for a in [dT_star_dt_mean, dT_star_dt_std, TOT_F_mean, TOT_F_std]
    a[mask] = missing_value
end

time = collect(Float64, 1:12)

filename = "data/dT_star_dt-TOT_F.nc"

NetCDFHelper.specialCopyNCFile(fn, filename, ["lat", "lon", "lat_vertices", "lon_vertices"])


for obj in [
    [
        T_star_mean, "T_star", Dict(
            "long_name"=>"Sea Surface Heat Content Density",
            "units"=>"J / m^2 / m",
            "missing_value" => missing_value
        )
    ], [
        dT_star_dt_mean, "dT_star_dt", Dict(
            "long_name"=>"Sea Surface Heat Content Density Changing Rate",
            "units"=>"W / m^2 / m",
            "missing_value" => missing_value
        )
    ], [
        dT_star_dt_std, "dT_star_dt_std", Dict(
            "long_name"=>"Sea Surface Heat Content Density Changing Rate Standard Deviation",
            "units"=>"W / m^2 / m",
            "missing_value" => missing_value
        )
    ], [
        TOT_F_mean, "TOT_F", Dict(
            "long_name"=>"Total Downward Traditional Heat Flux Mean",
            "units"=>"W / m^2",
            "missing_value" => missing_value
        )
    ], [
        TOT_F_std, "TOT_F_std", Dict(
            "long_name"=>"Total Downward Traditional Heat Flux Standard Deviation",
            "units"=>"W / m^2",
            "missing_value" => missing_value
        )
    ]
]

    var     = obj[1]
    varname = obj[2]
    varatts = obj[3]
    nccreate(
        filename,
        varname,
        "rlon",
        "rlat",
        "time", time,
        atts=varatts
    )
    ncwrite(var, filename, varname)

end

ncclose(filename)

@printf("Output file: %s\n", filename)
