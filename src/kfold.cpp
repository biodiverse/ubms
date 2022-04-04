// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::export]]
arma::mat get_loglik_occu(arma::vec y, int M, arma::imat si, arma::mat psimat,
                          arma::mat pmat, arma::ivec Kmin){

  int S = psimat.n_cols;
  mat out(S,M);

  vec psi_s(psimat.n_rows);
  vec p_s(pmat.n_rows);

  vec psite;
  double psi;
  int J;
  double cp;

  vec psub, ysub;

  for (unsigned s = 0; s < S; s++){

    psi_s = psimat.col(s);
    p_s = pmat.col(s);

    for (unsigned m = 0; m < M; m++){

      ysub = y.subvec(si(m,0), si(m,1));
      psi = psi_s(m);
      psub = p_s.subvec(si(m,0), si(m,1));
      J = psub.size();
      cp = 1.0;

      for (unsigned j = 0; j < J; j++){
        cp *= pow(psub(j), ysub(j)) * pow(1-psub(j), 1-ysub(j));
      }

      if(Kmin(m) == 1){
        out(s,m) = log(cp * psi + DBL_MIN);
      } else if(Kmin(m) == 0){
        out(s,m) = log(cp * psi + (1-psi) + DBL_MIN);
      }
    }
  }

  return out;
}
