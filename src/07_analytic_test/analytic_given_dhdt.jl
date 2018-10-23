using Printf
function prtArr(A)
    for i = 1:size(A)[1]
        for j = 1:size(A)[2]

            @printf("%.1f ", A[i,j]) 

        end
        @printf("\n")
    end
end


ρ    = 1027.0  # kg / m^3
c_p  = 3985.0  # J / kg / K

dtype = typeof(ρ)


years = 10
months_per_year = 12
pts_per_month = 1
pts_per_year = months_per_year * pts_per_month
secs_per_year = 86400.0 * 365.0

Ts_init = 300.00 
Td      = 0.0

θs_init = Ts_init * ρ * c_p
θd = Td * ρ * c_p

ret_len = months_per_year * pts_per_month * years

t = collect(range(0.0, stop=years * secs_per_year, length=ret_len+1))[1:end-1]
Δt = t[2] - t[1]

t_year = t / secs_per_year

S = t * 0.0
h = copy(S)
B = copy(S)

S0 = 100.0
for i = 1 : length(S)
    S[i] = S0 * sin(2π * t_year[i])
    h[i] = 30.0 #- 10.0 * sin(2π * t_year[i])
end

A = S0 * secs_per_year / (2π * 30.0)
c0 = A + θs_init
θs = - A * cos.(2π * t_year) .+ c0
Ts = θs / ρ / c_p


println(Ts)
