real lp_occu(array[] int y, real logit_psi, vector logit_p, int Kmin){
  real out;
  out = log_inv_logit(logit_psi) + bernoulli_logit_lpmf(y | logit_p);
  if(Kmin == 1){
    return out;
  }
  return log_sum_exp(out, log1m_inv_logit(logit_psi));
}

vector get_loglik_occu(array[] int y, int M, array[,] int J, array[,] int si, vector logit_psi,
                  vector logit_p, array[] int Kmin){
  vector[M] out;
  for (i in 1:M){
    out[i] = lp_occu(y[si[i,1]:si[i,2]], logit_psi[i],
                     logit_p[si[i,1]:si[i,2]], Kmin[i]);
   }
  return out;
}
