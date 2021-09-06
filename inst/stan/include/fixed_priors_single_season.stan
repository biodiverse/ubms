beta_state ~ normal(locations_state, scales_state);
beta_det ~ normal(locations_det, scales_det);
beta_scale ~ normal(locations_scale[1], scales_scale[1]);
beta_shape ~ normal(locations_shape[1], scales_shape[1]);
