int idx = 1;
if(has_random_state){
  for (i in 1:n_group_vars_state){
    b_state[idx:(n_random_state[i]+idx-1)] ~ normal(0, sigma_state[i]);
    idx += n_random_state[i];
  }
}

idx = 1;
if(has_random_det){
  for (i in 1:n_group_vars_det){
    b_det[idx:(n_random_det[i]+idx-1)] ~ normal(0, sigma_det[i]);
    idx += n_random_det[i];
  }
}
