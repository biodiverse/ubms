beta_state ~ normal(locations_state, scales_state);
beta_det ~ normal(locations_det, scales_det);
beta_scale ~ normal(0, 2.5);
beta_shape ~ normal(0, 2.5);
