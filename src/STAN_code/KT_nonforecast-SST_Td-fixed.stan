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
    real Q_std;
}

transformed data {
    int<lower=2*period+1> trimmed_N = raw_N - 2 * period;
    
    real true_future_theta[trimmed_N];

    // "_s" means staggered grid point
    real F_s[trimmed_N];
    real theta_s[trimmed_N];
    real theta[trimmed_N];

    if(raw_N % period != 0) {
       reject("raw_N = ", raw_N, " is not multiple of ", period);
    }

    // The following variables starts from index [period + 2] (e.g. February of 2nd year)
    // because we are making predictioins of next month. So the compared answer is shifted
    // by 1.
    true_future_theta = raw_theta[period+2:period+2+trimmed_N-1];
    theta             = raw_theta[period+1:period+1+trimmed_N-1];

    for (i in 1:trimmed_N) {
        F_s[i]     = (raw_F[i + period] + raw_F[i + period + 1]) / 2.0;
        theta_s[i] = (raw_theta[i + period] + raw_theta[i + period + 1]) / 2.0;
    }
}

parameters {
    real<lower=1, upper=5000> h[period];
    real Q_s[period];
}

model{

    real h_s[period];
    real dhdt_s[period];
    real we_s[period];

    real h_s_extended[trimmed_N];
    real we_s_extended[trimmed_N];
    real Q_s_extended[trimmed_N];
    real epsilon[trimmed_N];

    for(i in 1:period-1) {
        h_s[i] = (h[i] + h[i+1]) / 2.0;
    }
    h_s[period] = ( h[period] + h[1] ) / 2.0;


    for(i in 1:period-1) {
        dhdt_s[i] = (h[i+1] - h[i]) / dt;
    }
    dhdt_s[period] = ( h[1] - h[period] ) / dt;

    for(i in 1:period) { 
        we_s[i] = (dhdt_s[i] > 0) ? dhdt_s[i] : 0.0;
    }

    h_s_extended  = repeat_fill(h_s,  period, trimmed_N);
    we_s_extended = repeat_fill(we_s, period, trimmed_N);
    Q_s_extended  = repeat_fill(Q_s,  period, trimmed_N);

    for(i in 1:trimmed_N) {
        epsilon[i] = true_future_theta[i] - 
         ( theta[i] + (F_s[i] + Q_s_extended[i] - (theta_s[i] - theta_d) * we_s_extended[i]) / h_s_extended[i] * dt );
    }

    // Prior of Q
    for(i in 1:period) {
        Q_s[i] ~ normal(0, Q_std);
    }
    
    // Likelihood
    for(i in 1:trimmed_N) {
        epsilon[i] ~ normal(0, theta_std);
    }

}
