functions{

real lp_rn(int[] y, real log_lambda, vector logit_r, int J, int K, int Kmin){

  int numN = K - Kmin + 1;
  vector[J] q = 1 - inv_logit(logit_r);

  vector[numN] lp;
  vector[J] p;
  int N;

  for (i in 1:numN){
    N = K - i + 1;
    for (j in 1:J) p[j] = 1 - q[j]^N;
    lp[i] = poisson_log_lpmf(N | log_lambda) + bernoulli_lpmf(y | p);
  }
  return log_sum_exp(lp);
}

vector get_loglik_rn(int[] y, int M, int[,] J, int[,] si, vector log_lambda,
                     vector logit_p, int K, int[,] Kmin){
  vector[M] out;
  for (i in 1:M){
    out[i] = lp_rn(y[si[i,1]:si[i,2]], log_lambda[i], logit_p[si[i,1]:si[i,2]],
                   J[i,1], K, Kmin[i,1]);
  }
  return out;
}

}

data{

#include /include/data.stan

}

parameters{

#include /include/params_single_season.stan

}

transformed parameters{

vector[M] log_lambda;
vector[R] logit_p;
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

log_lik = get_loglik_rn(y, M, J, si, log_lambda, logit_p, K, Kmin);

}

model{

#include /include/rand_priors_single_season.stan
#include /include/fixed_priors_single_season.stan

target += sum(log_lik);

}

