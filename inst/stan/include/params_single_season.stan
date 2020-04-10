//Basic parameters for single-season model
vector[n_fixed_state] beta_state;
vector[n_fixed_det] beta_det;
vector<lower=0>[n_group_vars_state] sigma_state;
vector<lower=0>[n_group_vars_det] sigma_det;
vector[sum(n_random_state)] b_state;
vector[sum(n_random_det)] b_det;
