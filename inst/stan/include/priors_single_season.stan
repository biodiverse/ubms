// Fixed effects priors
target += lp_priors(beta_state, prior_dist_state, prior_pars_state);
target += lp_priors(beta_det, prior_dist_det, prior_pars_det);
target += lp_priors(beta_scale, prior_dist_scale, prior_pars_scale);
target += lp_priors(beta_shape, prior_dist_shape, prior_pars_shape);

// Random effects priors
target += lp_random_prior(has_random_state, n_group_vars_state, b_state,
                          n_random_state, sigma_state, prior_dist_state[3],
                          prior_pars_state);
target += lp_random_prior(has_random_det, n_group_vars_det, b_det,
                          n_random_det, sigma_det, prior_dist_det[3],
                          prior_pars_det);
