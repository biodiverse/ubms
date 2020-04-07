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
int occ_has_random;
int det_has_random;
int nFP_occ;
int nFP_det;
int n_grpvars_occ;
int n_grpvars_det;
int nRE_occ[occ_has_random ? n_grpvars_occ : 1];
int nRE_det[det_has_random ? n_grpvars_det: 1];
matrix[M, nFP_occ] X_occ;
matrix[M*J, nFP_det] X_det;
matrix[occ_has_random ? M : 0,sum(nRE_occ)] Z_occ;
matrix[det_has_random ? M*J : 0,sum(nRE_det)] Z_det;

}

parameters{

vector[nFP_occ] beta_occ;
vector[nFP_det] beta_det;
vector<lower=0>[n_grpvars_occ] sigma_occ;
vector<lower=0>[n_grpvars_det] sigma_det;
vector[sum(nRE_occ)] b_occ;
vector[sum(nRE_det)] b_det;

}

transformed parameters{

vector[M] logit_psi;
vector[M*J] logit_p;
vector[M] log_lik;

logit_psi = X_occ * beta_occ;
logit_p = X_det * beta_det;

if(occ_has_random){
  logit_psi = logit_psi + Z_occ * b_occ; 
}
if(det_has_random){
  logit_p = logit_p + Z_det * b_det;
}

log_lik = get_loglik(y, M, J, logit_psi, logit_p, no_detects);

}

model{

int idx = 1;

beta_occ ~ cauchy(0,2.5);
beta_det ~ cauchy(0,2.5);


if(occ_has_random){
  for (i in 1:n_grpvars_occ){
    b_occ[idx:(nRE_occ[i]+idx-1)] ~ normal(0, sigma_occ[i]);
    idx += nRE_occ[i];
  }
}

idx = 1;
if(det_has_random){
  for (i in 1:n_grpvars_det){
    b_det[idx:(nRE_det[i]+idx-1)] ~ normal(0, sigma_det[i]);
    idx += nRE_det[i];
  }
}

target += sum(log_lik);

}

