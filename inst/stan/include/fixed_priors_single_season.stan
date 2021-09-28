target += lp_priors(beta_state, prior_dist_state, prior_pars_state);
target += lp_priors(beta_det, prior_dist_det, prior_pars_det);
target += lp_priors(beta_scale, prior_dist_scale, prior_pars_scale);
target += lp_priors(beta_shape, prior_dist_shape, prior_pars_shape);
