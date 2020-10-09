vector ttd_prob_exp(vector y, vector log_lam, int[] delta){
  int J = num_elements(y);
  vector[J] e_lamt;
  real lam;
  for (j in 1:J){
    lam = exp(log_lam[j]);
    e_lamt[j] = pow(lam, delta[j]) * exp(-lam*y[j]);
  }
  return e_lamt;
}

vector ttd_prob_weib(vector y, vector log_lam, int[] delta, real log_k){
  int J = num_elements(y);
  vector[J] e_lamt;
  real k = exp(log_k);
  real lam;
  for (j in 1:J){
    lam = exp(log_lam[j]);
    e_lamt[j] = pow(k*lam*pow(lam*y[j], (k-1)), delta[j]) *
                exp(-1*pow(lam*y[j], k));
  }
  return e_lamt;
}

real lp_occuTTD(vector y, real logit_psi, vector log_lam,
                real log_k, int[] delta, int ydist){

  int J = num_elements(y);
  vector[J] e_lamt;
  real psi = inv_logit(logit_psi);
  real lik;

  if(ydist == 1){ //exponential
    e_lamt = ttd_prob_exp(y, log_lam, delta);
  } else if(ydist == 3){ //weibull
    e_lamt = ttd_prob_weib(y, log_lam, delta, log_k);
  }

  lik = psi * prod(e_lamt) + (1-psi) * (1-max(delta));
  return log(lik);
}

vector get_loglik_occuTTD(vector y, int M, int[,] si, vector logit_psi,
                          vector log_lam, real log_k, int[] delta, int ydist){
  vector[M] out;
  for (i in 1:M){
    out[i] = lp_occuTTD(y[si[i,1]:si[i,2]], logit_psi[i], log_lam[si[i,1]:si[i,2]],
                        log_k, delta[si[i,1]:si[i,2]], ydist);
  }
  return out;
}
