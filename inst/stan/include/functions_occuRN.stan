real lp_rn(int[] y, real log_lambda, vector logit_r, int J, int K, int Kmin){

  int numN = K - Kmin + 1;
  vector[J] q = 1 - inv_logit(logit_r);

  vector[numN] lp;
  vector[J] p;
  int N;

  for (i in 1:numN){
    N = K - i + 1;
    for (j in 1:J) p[j] = 1 - q[j]^N;
    lp[i] = poisson_log_lpmf(N | log_lambda) + bernoulli_lpmf(y | p);
  }
  return log_sum_exp(lp);
}

vector get_loglik_rn(int[] y, int M, int[,] J, int[,] si, vector log_lambda,
                     vector logit_p, int K, int[] Kmin){
  vector[M] out;
  for (i in 1:M){
    out[i] = lp_rn(y[si[i,1]:si[i,2]], log_lambda[i], logit_p[si[i,1]:si[i,2]],
                   J[i,1], K, Kmin[i]);
  }
  return out;
}

