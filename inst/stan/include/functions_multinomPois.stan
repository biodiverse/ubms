//Removal sampling
//p can be any length > 1
vector pi_removal(vector p){
  int J = num_elements(p);
  vector[J] pi_out;
  pi_out[1] = p[1];
  for (j in 2:J){
    pi_out[j] = pi_out[j-1] / p[j-1] * (1-p[j-1]) * p[j];
  }
  return pi_out;
}

//Double observer
//p must have 2 elements
vector pi_double(vector p){
  vector[3] pi_out;
  pi_out[1] = p[1] * (1 - p[2]);
  pi_out[2] = p[2] * (1 - p[1]);
  pi_out[3] = p[1] * p[2];
  return pi_out;
}

vector pi_fun(int pi_type, vector p, int J){
  vector[J] out;
  if(pi_type == 0){
    out = pi_double(p);
  } else if(pi_type == 1){
    out = pi_removal(p);
  } else {
    reject("Invalid pi function type");
  }
  return out;
}

real lp_multinomPois(int[] y, real log_lambda, vector logit_p, int pi_type){

  real loglik = 0.0;
  real lam = exp(log_lambda);
  int J = num_elements(y);
  int np = num_elements(logit_p);
  vector[np] p;
  vector[J] cp;

  for (i in 1:np){
    p[i] = inv_logit(logit_p[i]);
  }
  cp = pi_fun(pi_type, p, J);

  for (j in 1:J){
    loglik += poisson_lpmf(y[j] | lam * cp[j]);
  }
  return loglik;
}

vector get_loglik_multinomPois(int[] y, int M, int[,] si, vector log_lambda,
                               vector logit_p, int pi_type){

  vector[M] out;
  int J = num_elements(logit_p) / M; // use %/% in future version of stan
  int pstart = 1;
  int pend;
  for (i in 1:M){
    pend = pstart + J - 1;
    out[i] = lp_multinomPois(y[si[i,1]:si[i,2]], log_lambda[i],
                             logit_p[pstart:pend], pi_type);
    pstart += J;
  }
  return out;
}
