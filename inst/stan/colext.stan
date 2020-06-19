functions{

//can shortcut here I think
vector get_pY(int[] y, vector logit_p, int nd){
  vector[2] out;
  out[1] = nd;
  out[2] = exp(bernoulli_logit_lpmf(y | logit_p));
  return out;
}

matrix get_phi(vector phi_raw){
  return to_matrix(phi_raw, 2, 2, 0);
}

real lp_colext(int[] y, int T, int[] J, vector psi, matrix phi_raw, 
               vector p, int[] nd){

  int idx = 1;
  int end = idx + J[1] - 1;
  vector Dpt = get_pY(y[idx:end], p[idx:end], nd[1]);
  matrix phi_prod = diag_pre_multiply(Dpt, get_phi(phi_raw[1,]));
  idx += J[1];

  for (t in 2:(T-1)){
    end = idx + J[t] - 1;
    Dpt = get_pY(y[idx:end], p[idx:end], nd[t]);
    phi_prod *= diag_pre_multiply(Dpt, get_phi(phi_raw[t,]));
    idx += J[t];
  }

  end = idx + J[T] - 1;
  Dpt = get_pY(y[idx:end], p[idx:end], nd[T]);

  return log(dot_product(psi * phi _prod, Dpt));
}

//needs fixed
vector get_loglik_colext(int[] y, int M, int T, int[,] J, matrix psi_raw,
                  matrix phi_raw, vector logit_p){
  vector[M] out;
  int idx = 1;
  int end;
  for (i in 1:M){
    end = idx + J[1,i] - 1;
    out[i] = lp_occu(y[idx:end], logit_psi[i], logit_p[idx:end], nd[i]);
    idx += J[1,i];
   }
  return out;
}

}

data{

#include /include/data_single_season.stan

}

transformed data{

int no_detects[M];
for (m in 1:M){
  no_detects[m] = 1 - Kmin[m];
}

}

parameters{

#include /include/params_single_season.stan

}

transformed parameters{

vector[M] logit_psi;
vector[R] logit_p;
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

