functions{

#include /include/functions_occu.stan
#include /include/functions_occuRN.stan
#include /include/functions_pcount.stan
#include /include/functions_keyfuns.stan
#include /include/functions_distsamp.stan
#include /include/functions_multinomPois.stan
#include /include/functions_occuTTD.stan

}

data{

#include /include/data.stan

}

transformed data{

int ind_state;
int ind_det;
ind_state = prior_dist_state[1] == 0 ? 1 : 2;
ind_det = prior_dist_det[1] == 0 ? 1 : 2;

}

parameters{

#include /include/params_single_season.stan

}

transformed parameters{

vector[M] lp_state;
vector[n_obs_det] lp_det;
vector[M] log_lik;

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

#include /include/call_loglik_functions.stan

}

model{

#include /include/rand_priors_single_season.stan
#include /include/fixed_priors_single_season.stan

target += sum(log_lik);

}

