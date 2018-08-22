@printf("Running %s\n", basename(@__FILE__))

function prtArr(A)
    for i = 1:size(A)[1]
        for j = 1:size(A)[2]

            @printf("%d ", A[i,j]) 

        end
        @printf("\n")
    end
end

ρ    = 1027.0  # kg / m^3
c_p  = 3985.0  # J / kg / K
mon_secs = 365.0 / 12.0 * 86400.0
dt  = 1.0 * mon_secs

data_path = joinpath(dirname(@__FILE__), "..", "..", "data")
fn = joinpath(data_path, "SMART_Omon_GFDL-ESM2G_historical_r1i1p1_186101-200112.nc")


