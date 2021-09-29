functions{

#include /include/functions_occu.stan
#include /include/functions_occuRN.stan
#include /include/functions_pcount.stan
#include /include/functions_keyfuns.stan
#include /include/functions_distsamp.stan
#include /include/functions_multinomPois.stan
#include /include/functions_occuTTD.stan
#include /include/functions_priors.stan

}

data{

#include /include/data.stan

}

transformed data{

int include_scale;
int include_shape;

include_scale = prior_dist_scale[1] == 0 ? 0 : 1;
include_shape = prior_dist_shape[1] == 0 ? 0 : 1;
}

parameters{

#include /include/params_single_season.stan

}

transformed parameters{

vector[M] lp_state;
vector[n_obs_det] lp_det;
vector[M] log_lik;
real log_scale;
real log_shape;

lp_state = X_state * beta_state + offset_state;
lp_det = X_det * beta_det + offset_det;

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

log_scale = 0;
if(include_scale){
  log_scale = beta_scale[1];
}
log_shape = 0;
if(include_shape){
  log_shape = beta_shape[1];
}

#include /include/call_loglik_functions.stan

}

model{

#include /include/priors_single_season.stan

target += sum(log_lik);

}

