data {
    int<lower=0> N;
    real<lower=0> dt;
    real theta[N];
    real F[N];
}

transformed data {
    int yrs;
    if(yrs % 12 != 0) {
       reject("N is not multiple of 12"); 
    }
    yrs = N / 12;

    int trimmed_N = N - 24;
}


parameters {
    real<lower=0, upper=200> h[12];
    real Q_s[12];
    real Td;
}

transformed parameters {
    // suffix "_s" means staggered
    real dhdt_s[trimmed_N];
    real h_s[12];

    


    for(i in 1:11) {
        dhdt_s[i] = (h[i+1] - h[i] ) / dt;
    }
    dhdt_s[12] = ( h_s[1] - h_s[12] ) / dt;


}


model{
    for(i in 2:yrs-1) {

        // calculate residue
        Td = 



        if(y[i] == 1 || y[i] == 6) {
            target += log(p);
        } else {
            target += log((1.0 - 2.0 * p) / 4.0);
        }
    }
}
