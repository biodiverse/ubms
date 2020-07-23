functions{

#include /include/functions_keyfuns.stan
#include /include/functions_distsamp.stan

}

data{

#include /include/data.stan
int point;
vector[J[1,1] + 1] db;
int keyfun;
vector[R] conv_const;

}

transformed data{

 //Used in integration for hazard key function
real x_r[0];
int x_i[0];

}

parameters{

#include /include/params_single_season.stan

}

transformed parameters{

vector[M] lp_state;
vector[n_obs_det] lp_det;
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

log_lik = get_loglik_distsamp(y, M, db, si, lp_state, lp_det, z_dist,
                              beta_scale, point, keyfun, conv_const, x_r, x_i);

}

model{

#include /include/rand_priors_single_season.stan
#include /include/fixed_priors_single_season.stan

target += sum(log_lik);

}

