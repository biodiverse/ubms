functions{

real lp_occu(int[] y, real logit_psi, vector logit_p, int nd){
  real out;
  out = log_inv_logit(logit_psi) + bernoulli_logit_lpmf(y | logit_p);
  if(nd == 0){
    return out;
  }
  return log_sum_exp(out, log1m_inv_logit(logit_psi));
}

vector get_loglik_occu(int[,] y, int M, int J, vector logit_psi,
                  vector logit_p, int[] nd){
  vector[M] out;
  int idx = 1;
  for (i in 1:M){
    out[i] = lp_occu(y[i], logit_psi[i], logit_p[idx:(idx+J-1)],
             nd[i]);
    idx += J;
   }
  return out;
}

}

data{

#include /include/data_single_season.stan
int no_detects[M];

}

parameters{

#include /include/params_single_season.stan

}

transformed parameters{

vector[M] logit_psi;
vector[M*J] logit_p;
vector[M] log_lik;

logit_psi = X_state * beta_state;
logit_p = X_det * beta_det;

if(has_random_state){
  logit_psi = logit_psi + 
              csr_matrix_times_vector(Zdim_state[1], Zdim_state[2], Zw_state,
                                      Zv_state, Zu_state, b_state);
}
if(has_random_det){
  logit_p = logit_p + 
            csr_matrix_times_vector(Zdim_det[1], Zdim_det[2], Zw_det,
                                    Zv_det, Zu_det, b_det);
}

log_lik = get_loglik_occu(y, M, J, logit_psi, logit_p, no_detects);

}

model{

#include /include/model_single_season.stan

}

