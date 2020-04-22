// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace arma;

// [[Rcpp::export]]
arma::umat getZ_pcount(arma::umat y, arma::mat lam_post, arma::cube p_post, 
                      unsigned K, arma::uvec Kmin, arma::uvec kvals){

  int M = y.n_rows;
  int J = y.n_cols;
  int nsamples = lam_post.n_cols;

  vec kprob(K+1);
  double pp, bp;

  umat Zpost(M, nsamples);

  for (unsigned i=0; i < nsamples; i++){
    for (unsigned m=0; m < M; m++){
      kprob.zeros();
      for (unsigned k=Kmin(m); k < (K+1); k++){
        pp = R::dpois(k, lam_post(m, i), 0);
        bp = 1.0;
        for (unsigned j=0; j < J; j++){
          bp *= R::dbinom(y(m,j), k, p_post(m,j,i), 0);
        }
        kprob(k) = pp * bp; 
      }
      kprob = kprob / sum(kprob);
      Zpost(m,i) = Rcpp::RcppArmadillo::sample(kvals, 1, false, kprob)(0);
    }
  }
  
  return Zpost;

}
