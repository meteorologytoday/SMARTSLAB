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
    real raw_theta[raw_N];
    real raw_F[raw_N];
    real theta_std;
    real Q_std;
}

transformed data {
    int<lower=2*period+1> trimmed_N = raw_N - 2 * period;
    
    real true_future_theta[trimmed_N];
    real F[trimmed_N];
    real theta[trimmed_N];

    if(raw_N % period != 0) {
       reject("raw_N = ", raw_N, " is not multiple of ", period);
    }

    // The following variables starts from index [period + 2] (e.g. February of 2nd year)
    // because we are making predictioins of next month. So the compared answer is shifted
    // by 1.
    true_future_theta = raw_theta[period+2:period+2+trimmed_N-1];
    
    theta             = raw_theta[period+1:period+1+trimmed_N-1];
    F                 =     raw_F[period+1:period+1+trimmed_N-1];

}

parameters {
    real<lower=1, upper=5000> h[period];
    real Q[period];
    real<lower=1027.0*3985.0*263.15, upper=1027.0*3985.0*283.15> theta_d; // -10 ~ +10 deg C  
}

model{


    real h_extended[trimmed_N];
    real we_extended[trimmed_N];
    real Q_extended[trimmed_N];
    real epsilon[trimmed_N];

    real dhdt[period];
    real we[period];

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

    // Assume flat prior of theta_d : do nothing

    // Prior of Q
    for(i in 1:period) {
        Q[i] ~ normal(0, Q_std);
    }
    
    // Likelihood
    for(i in 1:trimmed_N) {
        epsilon[i] ~ normal(0, theta_std);
    }



}
