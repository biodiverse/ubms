real lp_single_prior(vector x, int dist, row_vector pars1,
                     row_vector pars2, row_vector pars3){
  real out = 0.0;
  if(dist == 1){
    out += normal_lpdf(x | pars1, pars2);
  } else if(dist == 2){
    out += uniform_lpdf(x | pars1, pars2);
  } else if(dist == 3){
    out += student_t_lpdf(x | pars1, pars2, pars3);
  } else if(dist == 4){
    out += logistic_lpdf(x | pars1, pars2);
  } else if(dist == 5){
    out += gamma_lpdf(x | pars1, pars2);
  } else if(dist == 6){
    out += double_exponential_lpdf(x | pars1, pars2);
  }
  return out;
}


real lp_priors(vector beta, int[] dist, matrix pars){

  int idx;
  real out = 0.0;
  int nb = num_elements(beta);
  if(nb == 0) return out;
  idx = dist[1] == 0 ? 1 : 2;

  // intercept
  out += lp_single_prior(beta[1:1], dist[1], pars[1,1:1],
                         pars[2,1:1], pars[3,1:1]);

  // regression coefficients
  out += lp_single_prior(beta[idx:nb], dist[2], pars[1,idx:nb],
                         pars[2,idx:nb], pars[3,idx:nb]);

  return out;
}

real lp_random_prior(int has_random, int n_group_vars, vector b,
                     int[] n_random, vector sigma, int dist, matrix pars){
  int idx = 1;
  real out = 0;
  int par_idx = cols(pars);
  row_vector[n_group_vars] rep_par1 = rep_row_vector(pars[1,par_idx], n_group_vars);
  row_vector[n_group_vars] rep_par2 = rep_row_vector(pars[2,par_idx], n_group_vars);
  row_vector[n_group_vars] rep_par3 = rep_row_vector(pars[3,par_idx], n_group_vars);
  if(has_random){
    out += lp_single_prior(sigma, dist, rep_par1, rep_par2, rep_par3);
    for (i in 1:n_group_vars){
      out += normal_lpdf(b[idx:(n_random[i]+idx-1)] | 0, sigma[i]);
      idx += n_random[i];
    }
  }
  return out;
}
