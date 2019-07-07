functions {
    real[] repeat_fill(real[] input, int period, int len) {
        real output[len];
        for (i in 1:len) {
            output[i] = input[(i-1) % period + 1];
        }

        return output;
    }

    real interpolate(real x, real x0, real y0, real x1, real y1) {
        return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
    }


    
}

data {
    int<lower=0> raw_N;     // raw data lengths (monthly)
    real dom[12];           // days of month
    int  steps[12];         // how many sub intervals from mid-ith-month to mid-(i+1)th-month
    real h_max;             // maximum depth it could get

    real raw_T[raw_N];      // Average temperature of each month
    real raw_F[raw_N];      // Average heat fluxes of each month

    real T_std;

    real c_sw;
    real rho_sw;
}

transformed data {
    int<lower=2*period+1> N = raw_N - 2 * 12;
    int years = N / 12;
    int sum_steps = sum(steps);
    int total_steps = years * sum_steps;

    real T[N];
    real true_future_T[N];

    real F[total_steps];

    real dt[12];
    real sub_dt[12];

    real heat_capacity_density = c_sw * rho_sw;

    if(raw_N % 12 != 0) {
       reject("raw_N = ", raw_N, " is not multiple of 12");
    }

    for(m in 1:12) {
        dt[m] = 86400.0 * (dom[1 + (m - 1) % 12] + dom[1 + (m+1 - 1) % 12]);
        sub_dt[m] = dt[m] / steps[m];
    }

    // The following variables starts from index [period + 2] (e.g. February of 2nd year)
    // because we are making predictioins of next month. So the compared answer is shifted
    // by 1.
    true_future_T = raw_T[12+2:12+2+N-1];
    T             = raw_T[12+1:12+1+N-1];

    // interpolate fluxes
    int t = 1;
    int t_month = 12+1;
    for(y in 1:years) {
        for(m in 1:12) {
            for(k in 1:steps[j]) {
                F[t] = interpolate(k-1.0, 0.0, raw_F[t_month], steps[m] - 1.0, raw_F[t_month+1]);
                t += 1;
            }
            t_month += 1;
        }
    }
}

parameters {
    real<lower=1> h[12];
    real Q[12];
    real<lower=-1.8> T_d; // degC  
}

model{

    real dhdt[12];
    real we[12];

    real _heat_capacity[sum_steps];
    real _Q[sum_steps];

    for(m in 1:12) {
        if(h[m] > h_max) {
            reject("Exceeds h_max at month ", m, ": ", h[m]);
        }
    }


    for(m in 1:11) {
        dhdt[m] = (h[m+1] - h[m]) / dt[m];
    }
    dhdt[12] = (h[1] - h[12]) / dt[12];

    for(m in 1:12) { 
        we[m] = (dhdt[m] > 0) ? dhdt[m] : 0.0;
    }


    int t_in_year = 1;
    for(m in 1:12) {
        for(k in 1:steps[j]) {
            if (m < 12) {
                _heat_capacity[t_in_year] = interpolate(k-1.0, 0.0, h[m], steps[m] - 1.0, h[m+1]) * heat_capacity_density;
                _Q[t_in_year] = interpolate(k-1.0, 0.0, Q[m], steps[m] - 1.0, Q[m+1]);
            } else {
                _heat_capacity[t_in_year] = interpolate(k-1.0, 0.0, h[12], steps[12] - 1.0, h[1]) * heat_capacity_density;
                _Q[t_in_year] = interpolate(k-1.0, 0.0, Q[12], steps[12] - 1.0, Q[1]);
            }
    
            t_in_year += 1;
        }
    }

    int t = 1;
    for(y in 1:years) {
        t_in_year = 1;
        for(m in 1:12) {
            real predict_T = T[t];
            for(k in 1:steps[i]) {

                predict_T += ( F[t] + _Q[t_in_year] - (T[t] - T_d) * we[m]) / _heat_capacity[t_in_year] * sub_dt[m];

                t_in_year += 1;
                t += 1;
            }

            // Likelihood
            (predict_T - true_future_T[t]) ~ normal(0, T_std); 
        }
    }

    
    // Assume flat prior of theta_d : do nothing

    /*
    // Prior of F + Q
    for(i in 1:period) {
        Q[i] ~ normal(0, Q_std);
    }
    */

}
