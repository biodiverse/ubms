real lp_pcount_pois(array[] int y, real log_lambda, vector logit_p, int K, int Kmin){

  real fac = 1;
  real ff = exp(log_lambda) * prod(1 - inv_logit(logit_p));
  real N;
  real ky;
  int numN = K - Kmin;
  for (i in 1:numN){
    N = K - i + 1;
    ky = 1;
    for (j in 1:size(y)){
      ky *= N / (N - y[j]);
    }
    fac = 1 + fac * ff * ky / N;
  }
  return  poisson_log_lpmf(Kmin | log_lambda) +
          binomial_logit_lpmf(y | Kmin, logit_p) +
          log(fac);
}

vector get_loglik_pcount(array[] int y, int M, array[,] int J, array[,] int si, vector log_lambda,
                         vector logit_p, int z_dist, real beta_scale, int K, array[] int Kmin){
  vector[M] out;
  for (i in 1:M){
    out[i] = lp_pcount_pois(y[si[i,1]:si[i,2]], log_lambda[i],
                            logit_p[si[i,1]:si[i,2]], K, Kmin[i]);
  }
  return out;
}

