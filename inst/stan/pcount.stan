functions{

real lp_pcount_pois(int[] y, real log_lambda, vector logit_p, int K, int Kmin){

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

vector get_loglik_pcount(int[,] y, int M, int J, vector log_lambda, vector logit_p, 
                         int mixture, real mix_param, int K, int[] Kmin){
  vector[M] out;
  int idx = 1;
  for (i in 1:M){
    out[i] = lp_pcount_pois(y[i], log_lambda[i], logit_p[idx:(idx+J-1)], K, Kmin[i]);
    idx += J;
  }
  return out;
}

}

data{

#include /include/data_single_season.stan
int K;
int Kmin[M];
int mixture;

}

parameters{
  
#include /include/params_single_season.stan
real beta_mix; //2nd parameter for abundance

}

transformed parameters{

vector[M] log_lambda;
vector[M*J] logit_p;
vector[M] log_lik;

log_lambda = X_state * beta_state;
logit_p = X_det * beta_det;

if(has_random_state){
  log_lambda = log_lambda + 
              csr_matrix_times_vector(Zdim_state[1], Zdim_state[2], Zw_state,
                                      Zv_state, Zu_state, b_state);
}
if(has_random_det){
  logit_p = logit_p + 
            csr_matrix_times_vector(Zdim_det[1], Zdim_det[2], Zw_det,
                                    Zv_det, Zu_det, b_det);
}

log_lik = get_loglik_pcount(y, M, J, log_lambda, logit_p, mixture, 
                            beta_mix, K, Kmin);

}

model{

#include /include/model_single_season.stan

}

