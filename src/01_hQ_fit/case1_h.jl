include("config.jl")

function nansum(A::Array)
    s, n = 0.0, 0.0
    for val in A
        if !isnan(val)
            s += val
            n += 1.0
        end
    end
    return s / n
end


# Discard the first and last year
TOT_F = TOT_F[:, :, 13:end-12] 
dT_star_dt = (T_star[:, :, 14:end-11] - T_star[:, :, 12:end-13]) / (2.0 * mon_secs)

#nyrs = 8
#TOT_F = TOT_F[:, :, 13:13+nyrs*12-1] 
#dT_star_dt = (T_star[:, :, 14:13+nyrs*12] - T_star[:, :, 12:13+nyrs*12-2]) / (2.0 * mon_secs)



print("Doing calculation... ")
h = nansum(TOT_F .* dT_star_dt) / nansum(dT_star_dt .^ 2.0)
println("done.")
@printf("The optimal mixed-layer thickness: %.2f m\n", h)

