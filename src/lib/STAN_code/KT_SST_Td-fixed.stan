functions {
    vector repeat_fill(vector input, int period, int len) {
        vector[len] output;
        for (i in 1:len) {
            output[i] = input[(i-1) % period + 1];
        }

        return output;
    }

    real dydx_parabolic(
        real x1, real y1,
        real x2, real y2,
        real x3, real y3
    ) {
    //
    //  This function returns the derivative at x_2 using
    //  parabolic function:
    // 
    //  (1) y = ax^2 + bx + c
    //
    //  from (1) we have its analytic derivative:
    //
    //      dy
    //  (2) --  =  2 ax + b
    //      dx
    //
    //  Suppose there are x_1, x_2, x_3 and
    //  corresponding y_1, y_2, y_3 then it
    //  can be easily proved that:
    //  
    //      dy|
    //  (3) --|    =  2 ax_2 + b
    //      dx|x_2
    //                y_1 - y_2     y_2 - y_3     y_1 - y_3
    //             =  ---------  +  ---------  -  ---------
    //                x_1 - x_2     x_2 - x_3     x_1 - x_3
    //

        return (y1 - y2)/(x1 - x2) + (y2 - y3)/(x2 - x3) - (y1 - y3)/(x1 - x3)
    }


    real KT_pred_theta(
        real theta,
        real theta_d,
        real F,
        real Q,
        real dhdt,
        real dt
    ) {
        return theta + (F + Q - (theta - theta_d) * ( (dhdt > 0) * dhdt ) ) * dt
    }

}

data {
    int<lower=0> raw_N;
    int<lower=1> period;
    real<lower=0> dt;
    real<lower=0> theta_d;
    vector[raw_N] raw_theta;
    vector[raw_N] raw_F;
}

transformed data {
    int<lower=2*period+1> trimmed_N = raw_N - 2 * period;
    
    vector[trimmed_N] true_future_theta;
    vector[trimmed_N] F;
    vector[trimmed_N] theta;

    real dtt = 2.0 * dt;
 
    if(all_N % period != 0) {
       reject("raw_N = ", raw_N, " is not multiple of ", period);
    }

    // The following variables starts from index [period + 2] (e.g. February of 2nd year)
    // because we are making predictioins of next month. So the compared answer is shifted
    // by 1.
    true_future_theta = raw_theta[period+2:period+2+trimmed_N-1]
    
    theta             = raw_theta[period+1:period+1+trimmed_N-1]
    F                 =     raw_F[period+1:period+1+trimmed_N-1]

}

parameters {
    real<lower=1, upper=5000> h[period];
    real Q[period];
}

model{


    vector[trimmed_N] h_extended;
    vector[trimmed_N] we_extended;
    vector[trimmed_N] Q_extended;
    vector[trimmed_N] epsilon;

    vector[period] dhdt;
    vector[period] we;

    for(i in 2:period-1) {
        dhdt[i] = (h[i+1] - h[i-1]) / dtt
    }
    dhdt[1]      = ( h[2] - h[period  ] ) / dtt;
    dhdt[period] = ( h[1] - h[period-1] ) / dtt;

    for(i in 1:period) { 
        we[i] = (dhdt[i] > 0) * dhdt[i]
    }

    h_extended  = repeat_fill(h,  period, trimmed_N);
    we_extended = repeat_fill(we, period, trimmed_N);
    Q_extended  = repeat_fill(Q,  period, trimmed_N);

    epsilon = true_future_theta -
        ( theta + (F + Q_extended - (theta - theta_d) * we) * dt )

    // Prior of Q
    for(i in 1:period) {
        Q[i] ~ normal(0, Q_std);
    }
    
    // Likelihood
    for(i in 1:trimmed_N) {
        epsilon[i] ~ normal(0, epsilon_std);
    }

}
