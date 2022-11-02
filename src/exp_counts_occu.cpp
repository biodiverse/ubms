// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;

double peh_occu(vec obs, int nd, vec p, double psi){
  int J = p.n_elem;
  double lik = psi;
  double cp;
  for (int j=0; j < J; j++){
    cp = p(j)*obs(j) + (1-p(j))*(1-obs(j));
    if(is_finite(cp)){
      lik *= cp;
    }
  }
  if(nd) lik += (1 - psi);
  return(lik);
}



// [[Rcpp::export]]
arma::vec exp_counts_occu(arma::mat obs, arma::ivec no_detects,
                             arma::vec psi, arma::vec p){

  int M = psi.n_elem;
  int J = p.n_elem / M;
  int n_eh = obs.n_cols;

  vec counts_expect = zeros(n_eh);
  int pstart, pend;
  for (int i=0; i < n_eh; i++){
    pstart = 0;
    for (int m=0; m < M; m++){
      pend = pstart + J - 1;
      counts_expect(i) += peh_occu(obs.col(i), no_detects(i),
                                   p.subvec(pstart, pend), psi(m));
      pstart += J;
    }
  }
  return(counts_expect);
}


double peh_occuRN(vec obs, int Kmin, vec r, double lam){
  int J = r.n_elem;
  double lik = 0.0;
  double obs_lik;
  double p;
  double cp;
  vec q = 1 - r;

  for (int k=Kmin; k<20; k++){
    obs_lik = 1.0;
    for (int j=0; j<J; j++){
      p = 1 - pow(q(j), k);
      cp = p*obs(j) + (1-p)*(1-obs(j));
      if(is_finite(cp)) obs_lik *= cp;
    }
    lik += R::dpois(k, lam, 0) * obs_lik;
  }

  return(lik);
}

// [[Rcpp::export]]
arma::vec exp_counts_occuRN(arma::mat obs, arma::ivec Kmin,
                             arma::vec lam, arma::vec r){

  int M = lam.n_elem;
  int J = r.n_elem / M;
  int n_eh = obs.n_cols;

  vec counts_expect = zeros(n_eh);
  int pstart, pend;
  for (int i=0; i < n_eh; i++){
    pstart = 0;
    for (int m=0; m < M; m++){
      pend = pstart + J - 1;
      counts_expect(i) += peh_occuRN(obs.col(i), Kmin(i),
                                   r.subvec(pstart, pend), lam(m));
      pstart += J;
    }
  }
  return(counts_expect);
}
