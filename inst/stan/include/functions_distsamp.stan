real lp_distsamp(int[] y, vector db, real log_lambda, real par1, real par2,
                 int point, int keyfun, vector conv_const){

  real lam = exp(log_lambda);
  int J = num_elements(db) - 1;
  real loglik = 0.0;
  real cp;
  for (j in 1:J){
    cp = prob_dist(par1, par2, keyfun, db[j], db[j+1], point);
    cp = cp * conv_const[j];
    loglik += poisson_lpmf(y[j] | lam * cp);
  }
  return loglik;
}

vector get_loglik_distsamp(int[] y, int M, vector db, int[,] si,
                           vector log_lambda, vector trans_par1, int z_dist,
                           real trans_par2, int point, int keyfun, vector conv_const){

  vector[M] out;
  for (i in 1:M){
    out[i] = lp_distsamp(y[si[i,1]:si[i,2]], db, log_lambda[i], trans_par1[i],
                         trans_par2, point, keyfun, conv_const[si[i,1]:si[i,2]]);
  }
  return out;
}
