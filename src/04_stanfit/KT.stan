functions {
    vector repeat_fill(vector input, int period, int len) {
        vector[len] output;
        for (i in 1:len) {
            output[i] = input[(i-1) % period + 1];
        }

        return output;
    }
}

data {
    int<lower=0> N;
    int<lower=1> period;
    real<lower=0> dt;
    vector[N] theta;
    vector[N] F;
    real<lower=0> epsilon_std;
    real<lower=0> Q_std;
}

transformed data {
    int yrs;
    int<lower=2*period+1> trimmed_N = N - 2 * period;
    vector[trimmed_N] F_s;
    vector[trimmed_N] theta_s;
    vector[trimmed_N] dthetadt_s;
 
    if(N % period != 0) {
       reject("N = ", N, " is not multiple of ", period);
    }
    yrs = N / period;

    for (i in 1:trimmed_N) {
        F_s[i]     = (F[i + period] + F[i + period + 1]) / 2.0;
        theta_s[i] = (theta[i + period] + theta[i + period + 1]) / 2.0;
    }
 
    for (i in 1:trimmed_N) {
        dthetadt_s[i] = (theta[i + period + 1] - theta[i + period]) / dt;
    }
    
}


parameters {
    real<lower=0, upper=200> h[period];
    real Q_s[period];
    real<lower=270*4000000, upper=300*4000000> theta_d;
}

model{

    // "_s" means staggered grid point

    vector[trimmed_N] we_s;     // entrainment velocity
    vector[trimmed_N] dhds_s;
    vector[trimmed_N] h_s;
    vector[trimmed_N] q_s;
    vector[trimmed_N] epsilon;

    for(i in 1:period-1) {
        h_s[i] = (h[i] + h[i+1] ) / 2.0;
    }
    h_s[period] = ( h[period] + h[1] ) / 2.0;

    // Calculate the change of MLD
    for(i in 1:period-1) {
        dhds_s[i] = (h[i+1] - h[i] ) / dt;
    }
    dhds_s[period] = ( h[1] - h[period] ) / dt;

    // Entrainment condition
    for(i in 1:period) {
        we_s[i] = ( dhds_s[i] > 0 ) ? dhds_s[i] : 0;
        q_s[i] = Q_s[i];
    }

    we_s   = repeat_fill(we_s, period, trimmed_N);
    h_s    = repeat_fill(h_s, period, trimmed_N);
    q_s    = repeat_fill(q_s, period, trimmed_N);

    epsilon = h_s .* dthetadt_s + we_s .* (theta_s - theta_d) - F_s - q_s;

    // Prior of Q
    for(i in 1:period) {
        Q_s[i] ~ normal(0, Q_std);
    }

    // Likelihood
    for(i in 1:trimmed_N) {
        epsilon[i] ~ normal(0, epsilon_std);
    }

}
