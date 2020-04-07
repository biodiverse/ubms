functions{

real lp_occu(int[] y, real logit_psi, vector logit_p, int nd){
  real out;
  out = log_inv_logit(logit_psi) + bernoulli_logit_lpmf(y | logit_p);
  if(nd == 0){
    return out;
  }
  return log_sum_exp(out, log1m_inv_logit(logit_psi));
}
vector get_loglik(int[,] y, int M, int J, vector logit_psi,
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

int M;
int J;
int y[M,J];
int no_detects[M];
int has_random_state;
int has_random_det;
int n_fixed_state;
int n_fixed_det;
int n_group_vars_state;
int n_group_vars_det;
int n_random_state[has_random_state ? n_group_vars_state : 1];
int n_random_det[has_random_det ? n_group_vars_det: 1];
matrix[M, n_fixed_state] X_state;
matrix[M*J, n_fixed_det] X_det;
matrix[has_random_state ? M : 0,sum(n_random_state)] Z_state;
matrix[has_random_det ? M*J : 0,sum(n_random_det)] Z_det;

}

parameters{

vector[n_fixed_state] beta_state;
vector[n_fixed_det] beta_det;
vector<lower=0>[n_group_vars_state] sigma_state;
vector<lower=0>[n_group_vars_det] sigma_det;
vector[sum(n_random_state)] b_state;
vector[sum(n_random_det)] b_det;

}

transformed parameters{

vector[M] logit_psi;
vector[M*J] logit_p;
vector[M] log_lik;

logit_psi = X_state * beta_state;
logit_p = X_det * beta_det;

if(has_random_state){
  logit_psi = logit_psi + Z_state * b_state; 
}
if(has_random_det){
  logit_p = logit_p + Z_det * b_det;
}

log_lik = get_loglik(y, M, J, logit_psi, logit_p, no_detects);

}

model{

int idx = 1;

beta_state ~ cauchy(0,2.5);
beta_det ~ cauchy(0,2.5);

if(has_random_state){
  for (i in 1:n_group_vars_state){
    b_state[idx:(n_random_state[i]+idx-1)] ~ normal(0, sigma_state[i]);
    idx += n_random_state[i];
  }
}

idx = 1;
if(has_random_det){
  for (i in 1:n_group_vars_det){
    b_det[idx:(n_random_det[i]+idx-1)] ~ normal(0, sigma_det[i]);
    idx += n_random_det[i];
  }
}

target += sum(log_lik);

}

