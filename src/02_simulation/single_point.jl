
function f(
    T_star    :: G,
    F    :: G,
    Q    :: G,  
    h    :: G,
    dhdt :: G
) where G <: AbstractFloat

    @printf("T_star: %.2f, F: %.2f, Q: %.2f, h: %.2f, dhdt: %.2e\n", T_star, F, Q, h, dhdt)

    return ((F + Q) - T_star * dhdt) / h

end


function extend(a::AbstractArray, len::Int)
    if length(a) < len
        a = repeat(a, outer=ceil(Int, len / length(a)))
    end

    return a[1:len]
end


function linear_interpolate(a::Array{G, 1}, i::G) where G <: AbstractFloat
    return (i - floor(i)) * a[ceil(Int, i)] + (ceil(i) - i) * a[floor(Int, i)]
end


function simulate_single_point(
    init_T_star  :: G,
    dt      :: G,
    h       :: Array{G, 1},
    dhdt    :: Array{G, 1},
    Q       :: Array{G, 1},
    F       :: Array{G, 1},
    h_min   :: G,
    ret_len :: Int
) where G <: AbstractFloat


    T_star = zeros(G, ret_len)
    T_star[1] = init_T_star

    Q    = extend(Q,    ret_len)
    F    = extend(F,    ret_len)
    h    = extend(h,    ret_len)
    dhdt = extend(dhdt, ret_len)

    println(dhdt)

    for i=1:ret_len-1

        #=  RK4
        k1 = dt * f(t, y)
        k2 = dt * f(t+dt/2, y+k1/2)
        k3 = dt * f(t+dt/2, y+k2/2)
        k4 = dt * f(t+dt, y+k3)
        =#

        k1 = dt * f(
            T_star[i],
            F[i],
            Q[i],
            h[i],
            dhdt[i]
        )

        k2 = dt * f(
            T_star[i] + k1 / 2.0,
            (F[i]    + F[i+1])/2.0,
            (Q[i]    + Q[i+1])/2.0,
            (h[i]    + h[i+1])/2.0,
            (dhdt[i] + dhdt[i+1])/2.0
        )

        k3 = dt * f(
            T_star[i] + k2 / 2.0,
            (F[i]    + F[i+1])/2.0,
            (Q[i]    + Q[i+1])/2.0,
            (h[i]    + h[i+1])/2.0,
            (dhdt[i] + dhdt[i+1])/2.0
        )

        k4 = dt * f(
            T_star[i] + k3,
            F[i+1],
            Q[i+1],
            h[i+1],
            dhdt[i+1]
        )

        @printf("[Step %04d] %.2f, %.2f, %.2f, %.2f\n", i, k1, k2, k3, k4)

        T_star[i+1] = T_star[i] + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    end 

    return T_star
end

#=
        k2 = dt * f(
            T[i] + k1 / 2.0,
            (F[i] + F[i+1])/2.0,
            (F[i] + F[i+1])/2.0,
            (F[i] + F[i+1])/2.0,
            (F[i] + F[i+1])/2.0,
            linear_interpolate(F,    i+0.5),
            linear_interpolate(Q,    i+0.5),
            linear_interpolate(h,    i+0.5),
            linear_interpolate(dhdt, i+0.5)
        )

        k3 = dt * f(
            T[i] + k2 / 2.0,
            linear_interpolate(F,    i+0.5),
            linear_interpolate(Q,    i+0.5),
            linear_interpolate(h,    i+0.5),
            linear_interpolate(dhdt, i+0.5)
        )

        k4 = dt * f(
            T[i] + k3,
            linear_interpolate(F,    i+1.0),
            linear_interpolate(Q,    i+1.0),
            linear_interpolate(h,    i+1.0),
            linear_interpolate(dhdt, i+1.0)
        )

=#
