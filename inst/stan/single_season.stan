functions{

#include /include/functions_occu.stan
#include /include/functions_occuRN.stan
#include /include/functions_pcount.stan

}

data{

#include /include/data.stan

}

parameters{

#include /include/params_single_season.stan

}

transformed parameters{

vector[M] lp_state;
vector[R] lp_det;
vector[M] log_lik;

lp_state = X_state * beta_state;
lp_det = X_det * beta_det;

if(has_random_state){
  lp_state = lp_state +
              csr_matrix_times_vector(Zdim_state[1], Zdim_state[2], Zw_state,
                                      Zv_state, Zu_state, b_state);
}
if(has_random_det){
  lp_det = lp_det +
            csr_matrix_times_vector(Zdim_det[1], Zdim_det[2], Zw_det,
                                    Zv_det, Zu_det, b_det);
}

if(model_code == 0){
  log_lik = get_loglik_occu(y, M, J, si, lp_state, lp_det, Kmin[,1]);
} else if(model_code == 1){
  log_lik = get_loglik_rn(y, M, J, si, lp_state, lp_det, K, Kmin[,1]);
} else if(model_code == 2){
  log_lik = get_loglik_pcount(y, M, J, si, lp_state, lp_det, z_dist,
                              beta_scale, K, Kmin[,1]);
}

}

model{

#include /include/rand_priors_single_season.stan
#include /include/fixed_priors_single_season.stan

target += sum(log_lik);

}

