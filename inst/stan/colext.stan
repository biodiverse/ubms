functions{

#include /include/functions_priors.stan

//can shortcut here I think
vector get_pY(int[] y, vector logit_p, int nd){
  vector[2] out;
  out[1] = nd;
  out[2] = exp(bernoulli_logit_lpmf(y | logit_p));
  return out;
}

matrix phi_matrix(row_vector phi_raw){
  return to_matrix(phi_raw, 2, 2, 0);
}

//delta-step transition prob matrix via Chapman-Kolmogorov equation
matrix get_phi(matrix phi_raw, int Tstart, int Tnext){
  matrix[2,2] phi = diag_matrix(rep_vector(1, 2));
  int delta = Tnext - Tstart;
  if(delta == 1){
    return phi_matrix(phi_raw[Tstart,]);
  }

  for (d in 1:delta){
    phi = phi * phi_matrix(phi_raw[(Tstart + d - 1),]);
  }
  return phi;
}

//Ts = indices of primary periods when site was sampled (eg not all NA)
real lp_colext(int[] y, int[] Tsamp, int[] J, row_vector psi, matrix phi_raw,
               vector logit_p, int[] nd){

  int T = size(Tsamp);
  matrix[2,2] phi_prod = diag_matrix(rep_vector(1, 2));
  matrix[2,2] phi;
  vector[2] Dpt;
  int idx = 1;
  int end;

  if(T > 1){
    for (t in 1:(T-1)){
      phi = get_phi(phi_raw, Tsamp[t], Tsamp[t+1]);
      end = idx + J[t] - 1;
      Dpt = get_pY(y[idx:end], logit_p[idx:end], nd[t]);
      phi_prod *= diag_pre_multiply(Dpt, phi);
      idx += J[t];
    }
  }

  end = idx + J[T] - 1;
  Dpt = get_pY(y[idx:end], logit_p[idx:end], nd[T]);

  return log(dot_product(psi * phi_prod, Dpt));
}

//needs fixed
vector get_loglik_colext(int[] y, int M, int[] Tsamp, int[,] J, int[,] si,
                         matrix psi_raw, matrix phi_raw, vector logit_p,
                         int[,] nd){
  vector[M] out;
  for (i in 1:M){
    out[i] = lp_colext(y[si[i,1]:si[i,2]], Tsamp[si[i,3]:si[i,4]], J[i,],
                       psi_raw[i,], phi_raw[si[i,5]:si[i,6],],
                       logit_p[si[i,1]:si[i,2]], nd[i,]);
  }
  return out;
}

}

data{

#include /include/data.stan

int has_random_col;
int has_random_ext;
int n_obs_col;
int n_obs_ext;
int n_fixed_col;
int n_fixed_ext;
int n_group_vars_col;
int n_group_vars_ext;
int n_random_col[has_random_col ? n_group_vars_col : 1];
int n_random_ext[has_random_ext ? n_group_vars_ext: 1];
matrix[M*(T-1), n_fixed_col] X_col;
matrix[M*(T-1), n_fixed_ext] X_ext;
vector[M*(T-1)] offset_col;
vector[M*(T-1)] offset_ext;

int Zdim_col[5];
vector[Zdim_col[3]] Zw_col;
int Zv_col[Zdim_col[4]];
int Zu_col[Zdim_col[5]];

int Zdim_ext[5];
vector[Zdim_ext[3]] Zw_ext;
int Zv_ext[Zdim_ext[4]];
int Zu_ext[Zdim_ext[5]];

int prior_dist_col[3];
int prior_dist_ext[3];
matrix[3, (n_fixed_col+1)] prior_pars_col;
matrix[3, (n_fixed_ext+1)] prior_pars_ext;
}

transformed data{

int no_detects[M, T];
int include_scale;
int include_shape;

for (m in 1:M){
  for (t in 1:T){
    no_detects[m, t] = 1 - Kmin[m, t];
  }
}

include_scale = 0;
include_shape = 0;
}

parameters{

#include /include/params_single_season.stan

vector[n_fixed_col] beta_col;
vector[n_fixed_ext] beta_ext;
vector<lower=0>[n_group_vars_col] sigma_col;
vector<lower=0>[n_group_vars_ext] sigma_ext;
vector[sum(n_random_col)] b_col;
vector[sum(n_random_ext)] b_ext;

}

transformed parameters{

vector[M] logit_psi;
matrix[M,2] psi_raw;
vector[M*(T-1)] logit_col;
vector[M*(T-1)] logit_ext;
matrix[M*(T-1), 4] phi_raw;
vector[R] logit_p;
vector[M] log_lik;

//psi
logit_psi = X_state * beta_state + offset_state;
if(has_random_state){
  logit_psi = logit_psi +
              csr_matrix_times_vector(Zdim_state[1], Zdim_state[2], Zw_state,
                                      Zv_state, Zu_state, b_state);
}

for (i in 1:M){
  psi_raw[i,2] = inv_logit(logit_psi[i]);
  psi_raw[i,1] = 1 - psi_raw[i,2];
}

//phi
logit_col = X_col * beta_col + offset_col;
if(has_random_col){
  logit_col = logit_col +
              csr_matrix_times_vector(Zdim_col[1], Zdim_col[2], Zw_col,
                                      Zv_col, Zu_col, b_col);
}

logit_ext = X_ext * beta_ext + offset_ext;
if(has_random_ext){
  logit_ext = logit_ext +
              csr_matrix_times_vector(Zdim_ext[1], Zdim_ext[2], Zw_ext,
                                      Zv_ext, Zu_ext, b_ext);
}

for (i in 1:(M*(T-1))){
  phi_raw[i,2] = inv_logit(logit_col[i]);
  phi_raw[i,1] = 1 - phi_raw[i,2];
  phi_raw[i,3] = inv_logit(logit_ext[i]);
  phi_raw[i,4] = 1 - phi_raw[i,3];
}

//det
logit_p = X_det * beta_det + offset_det;
if(has_random_det){
  logit_p = logit_p +
            csr_matrix_times_vector(Zdim_det[1], Zdim_det[2], Zw_det,
                                    Zv_det, Zu_det, b_det);
}

log_lik = get_loglik_colext(y, M, Tsamp, J, si, psi_raw, phi_raw,
                            logit_p, no_detects);

}

model{

#include /include/priors_single_season.stan

target += lp_priors(beta_col, prior_dist_col, prior_pars_col);
target += lp_priors(beta_ext, prior_dist_ext, prior_pars_ext);

target += lp_random_prior(has_random_col, n_group_vars_col, b_col,
                          n_random_col, sigma_col, prior_dist_col[3],
                          prior_pars_col);
target += lp_random_prior(has_random_ext, n_group_vars_ext, b_ext,
                          n_random_ext, sigma_ext, prior_dist_ext[3],
                          prior_pars_ext);

target += sum(log_lik);

}

