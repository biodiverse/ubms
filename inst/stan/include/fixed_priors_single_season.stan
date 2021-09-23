if(prior_dist_state[1] == 1){
  beta_state[1] ~ normal(prior_pars_state[1,1], prior_pars_state[2,1]);
} else if(prior_dist_state[1] == 2){
  beta_state[1] ~ uniform(prior_pars_state[1,1], prior_pars_state[2,1]);
} else if(prior_dist_state[1] == 3){
  beta_state[1] ~ student_t(prior_pars_state[3,1], prior_pars_state[1,1],
                            prior_pars_state[2,1]);
}

if(prior_dist_state[2] == 1){
  beta_state[ind_state:n_fixed_state] ~ normal(prior_pars_state[1,ind_state:n_fixed_state],
                                               prior_pars_state[2,ind_state:n_fixed_state]);
} else if(prior_dist_state[2] == 2){
  beta_state[ind_state:n_fixed_state] ~ uniform(prior_pars_state[1,ind_state:n_fixed_state],
                                                prior_pars_state[2,ind_state:n_fixed_state]);
} else if(prior_dist_state[2] == 3){
  beta_state[ind_state:n_fixed_state] ~ student_t(prior_pars_state[3,ind_state:n_fixed_state],
                                                  prior_pars_state[1,ind_state:n_fixed_state],
                                                  prior_pars_state[2,ind_state:n_fixed_state]);
}

// Detection priors
if(prior_dist_det[1] == 1){
  beta_det[1] ~ normal(prior_pars_det[1,1], prior_pars_det[2,1]);
} else if(prior_dist_det[1] == 2){
  beta_det[1] ~ uniform(prior_pars_det[1,1], prior_pars_det[2,1]);
} else if(prior_dist_det[1] == 3){
  beta_det[1] ~ student_t(prior_pars_det[3,1], prior_pars_det[1,1], prior_pars_det[2,1]);
}

if(prior_dist_det[2] == 1){
  beta_det[ind_det:n_fixed_det] ~ normal(prior_pars_det[1,ind_det:n_fixed_det],
                                         prior_pars_det[2,ind_det:n_fixed_det]);
} else if(prior_dist_det[2] == 2){
  beta_det[ind_det:n_fixed_det] ~ uniform(prior_pars_det[1,ind_det:n_fixed_det],
                                          prior_pars_det[2,ind_det:n_fixed_det]);
} else if(prior_dist_det[2] == 3){
  beta_det[ind_det:n_fixed_det] ~ student_t(prior_pars_det[3,ind_det:n_fixed_det],
                                            prior_pars_det[1,ind_det:n_fixed_det],
                                            prior_pars_det[2,ind_det:n_fixed_det]);
}


if(prior_dist_scale[1] == 0){
  beta_scale ~ normal(0, 0.01);
} else if(prior_dist_scale[1] == 1){
  beta_scale ~ normal(prior_pars_scale[1,1], prior_pars_scale[2,1]);
} else if(prior_dist_scale[1] == 2){
  beta_scale ~ uniform(prior_pars_scale[1,1], prior_pars_scale[2,1]);
} else if(prior_dist_scale[1] == 3){
  beta_scale ~ student_t(prior_pars_scale[3,1], prior_pars_scale[1,1], prior_pars_scale[2,1]);
}

if(prior_dist_shape[1] == 0){
  beta_shape ~ normal(0, 0.01);
} else if(prior_dist_shape[1] == 1){
  beta_shape ~ normal(prior_pars_shape[1,1], prior_pars_shape[2,1]);
} else if(prior_dist_shape[1] == 2){
  beta_shape ~ uniform(prior_pars_shape[1,1], prior_pars_shape[2,1]);
} else if(prior_dist_shape[1] == 3){
  beta_shape ~ student_t(prior_pars_shape[3,1], prior_pars_shape[1,1], prior_pars_shape[2,1]);
}
