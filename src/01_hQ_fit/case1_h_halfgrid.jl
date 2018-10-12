include("config.jl")

@inline nansum(A::Array) = sum(A[isfinite.(A)])

# Discard the first and last year
TOT_F = (TOT_F[:, :, 13:end-12] + TOT_F[:, :, 14:end-11]) / 2.0 
dT_star_dt = (T_star[:, :, 14:end-11] - T_star[:, :, 13:end-12]) / mon_secs

#nyrs = 8
#TOT_F = TOT_F[:, :, 13:13+nyrs*12-1] 
#dT_star_dt = (T_star[:, :, 14:13+nyrs*12] - T_star[:, :, 12:13+nyrs*12-2]) / (2.0 * mon_secs)



print("Doing calculation... ")
h = nansum(TOT_F .* dT_star_dt) / nansum(dT_star_dt .^ 2.0)
println("done.")
@printf("The optimal mixed-layer thickness: %.2f m\n", h)


# Doing measure

ϵ = h * dT_star_dt - TOT_F
ϵ2sum = nansum(ϵ.^2)
count =  sum(isfinite.(ϵ))

@printf("Total Residue: %.e\n", ϵ2sum)
@printf("Valid count: %d\n", count)
@printf("Mean residue: %.e\n", (ϵ2sum / count)^0.5)


