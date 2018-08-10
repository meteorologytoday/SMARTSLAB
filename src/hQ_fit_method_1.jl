include("config.jl")
include("NetCDFHelper.jl")

import GLM
import DataFrames
using NetCDF



# Discard the first and last year
TOT_F = TOT_F[:, :, 13:end-12] 
dT_star = (SST[:, :, 14:end-11] - SST[:, :, 12:end-13]) * ρ * c_p


#N = size(TOT_F)[3]/12
#ϕ = zeros(eltype(TOT_F), N, 2)
#ϕ[:, 2] = 1.0

h = zeros(eltype(TOT_F), length(rlons), length(rlats), 12)
Q = copy(h)

h_std = copy(h)
Q_std = copy(h)


for m = 1:12
    println("Doing month [$m]")
    for i = 1:length(rlons), j = 1:length(rlats)

        if mask[i,j]
            h[i,j,m] = NaN
            Q[i,j,m] = NaN
            h_std[i,j,m] = NaN
            Q_std[i,j,m] = NaN
            continue
        end

        data = DataFrames.DataFrame(X=TOT_F[i,j,m:12:end], Y=dT_star[i,j,m:12:end])
        ols = GLM.lm(GLM.@formula(Y ~ X), data)

        a = GLM.coef(ols)[2]
        b = GLM.coef(ols)[1]

        a_std  = GLM.stderror(ols)[2]
        b_std  = GLM.stderror(ols)[1]
        ab_cov = GLM.vcov(ols)[1,2]

        h[i,j,m] = 2.0 * dt / a
        Q[i,j,m] = b / a

        h_std[i,j,m] = abs(h[i,j,m] * a_std / a)
        Q_std[i,j,m] = abs(Q[i,j,m]) * ((a_std / a)^2.0 + (b_std / b)^2.0 - 2.0 * ab_cov / (a * b))^(0.5)

    end
end

mask = isnan.(h)

h[mask] = missing_value
Q[mask] = missing_value

h_std[mask] = missing_value
Q_std[mask] = missing_value



time = collect(Float64, 1:12)

filename = "hQ_method_1.nc"
NetCDFHelper.specialCopyNCFile(fn, filename, ["lat", "lon", "lat_vertices", "lon_vertices"])


for obj in [
    [
        h, "h", Dict(
            "long_name"=>"Constant Mixed-layer Thickness",
            "units"=>"m",
            "missing_value" => missing_value
        )
    ], [
        Q, "Q", Dict(
            "long_name"=>"Q-flux",
            "units"=>"J / m^2",
            "missing_value" => missing_value
        )
    ], [
        h_std, "h_std", Dict(
            "long_name"=>"Constant Mixed-layer Thickness Standard Deviation",
            "units"=>"m",
            "missing_value" => missing_value
        )
    ], [
        Q_std, "Q_std", Dict(
            "long_name"=>"Q-flux Standard Deviation",
            "units"=>"J / m^2",
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
