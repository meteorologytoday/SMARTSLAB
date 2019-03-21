functions {
    real[] repeat_fill(real[] input, int period, int len) {
        real output[len];
        for (i in 1:len) {
            output[i] = input[(i-1) % period + 1];
        }

        return output;
    }
}

data {

    int<lower=0> raw_N;
    int<lower=1> period;
    real<lower=0> dt;
    real<lower=0> theta_d;
    real raw_theta[raw_N];
    real raw_F[raw_N];
    real theta_std;
    real theta_trend_std;
}

transformed data {
    int<lower=2*period+1> trimmed_N = raw_N - 2 * period;
    
    real true_future_theta[trimmed_N];
    real F[trimmed_N];
    real theta[trimmed_N];

    real F_trend;
    int period_N;


    if(raw_N % period != 0) {
       reject("raw_N = ", raw_N, " is not multiple of ", period);
    }


    period_N = raw_N / period;
    // The following variables starts from index [period + 2] (e.g. February of 2nd year)
    // because we are making predictioins of next month. So the compared answer is shifted
    // by 1.
    true_future_theta = raw_theta[period+2:period+2+trimmed_N-1];
    
    theta             = raw_theta[period+1:period+1+trimmed_N-1];
    F                 =     raw_F[period+1:period+1+trimmed_N-1];

    F_trend = sum(F) / (period_N * period * dt);
}

parameters {
    real<lower=1, upper=2000> h[period];
    real Q[period];
}

model{


    real h_extended[trimmed_N];
    real we_extended[trimmed_N];
    real Q_extended[trimmed_N];
    real epsilon[trimmed_N];

    real dhdt[period];
    real we[period];

    real total_trend;
    real mean_h;

    total_trend = sum(Q) / (dt * period) + F_trend;
    mean_h = sum(h) / period;

    for(i in 2:period) {
        dhdt[i] = (h[i] - h[i-1]) / dt;
    }
    dhdt[1] = ( h[1] - h[period  ] ) / dt;

    for(i in 1:period) { 
        we[i] = (dhdt[i] > 0) ? dhdt[i] : 0.0;
    }

    h_extended  = repeat_fill(h,  period, trimmed_N);
    we_extended = repeat_fill(we, period, trimmed_N);
    Q_extended  = repeat_fill(Q,  period, trimmed_N);

    for(i in 1:trimmed_N) {
        epsilon[i] = true_future_theta[i] -
         ( theta[i] + (F[i] + Q_extended[i] - (theta[i] - theta_d) * we_extended[i]) / h_extended[i] * dt );
    }

    // Prior of climate drift
    total_trend ~ normal(0, theta_trend_std);

    /*
    // Prior of Q
    for(i in 1:period) {
    }
    */
    
    // Likelihood
    for(i in 1:trimmed_N) {
        epsilon[i] ~ normal(0, theta_std);
    }

}
