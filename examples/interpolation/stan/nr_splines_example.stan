functions {
  #include nr_splines.stanfunctions
}
data {

  int<lower=2> N_in;		// number of incoming datapoints
  int<lower=1> N_out;		//   and outgoing

  vector[N_in] x;		// knots
  vector[N_in] y;		//   and matching values
  vector[N_out] x_out;		// interpolant points
  vector[N_out] y_data;		//   and noisy data points to fit

}
parameters {

  vector[N_in] y_est;		// estimated original values
  real<lower=0> stdev;          // standard deviation

}
transformed parameters {

  // calculate second derivatives, with no constraint on endpoints
  vector[N_in] secder = nr__second_derivative( x, y_est, 1e32, 1e32 );

  // interpolate
  vector[N_out] y_hat = nr__interpolate( x_out, x, y_est, secder );

  // gotta do derivatives manually
  vector[N_out] yp_hat;
  { // trick Stan into allowing integer variables in this block
  int right = 2;	

  for( j in 1:N_out ) {

	while( (x[right] < x_out[j]) && (right < N_in) ) 
		right += 1; 

	yp_hat[j] = nr__derivative( x_out[j], x[right-1], x[right], y_est[right-1], y_est[right], secder[right-1], secder[right] );

	}
  }
}
model {

  // priors
  y_est ~ normal( 0, 2 );
  stdev ~ exponential( .1 );

  // likelihood
  y_hat ~ normal( y_data, stdev );

}
generated quantities {

  // assess the quality of the fit
  real ssquares = sum( (y_data - y_hat)^2 );

}
