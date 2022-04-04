// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::export]]
arma::mat get_loglik_occu(arma::vec y, int M, arma::imat si, arma::mat psimat,
                          arma::mat pmat, arma::ivec Kmin){

  int S = psimat.n_cols;
  mat out(S,M);
  vec p_s(pmat.n_rows);
  int J;
  double cp;
  vec psub, ysub;

  for (int s = 0; s < S; s++){

    p_s = pmat.col(s);

    for (int m = 0; m < M; m++){

      ysub = y.subvec(si(m,0), si(m,1));
      psub = p_s.subvec(si(m,0), si(m,1));
      J = psub.size();
      cp = 1.0;

      for (int j = 0; j < J; j++){
        cp *= pow(psub(j), ysub(j)) * pow(1-psub(j), 1-ysub(j));
      }

      if(Kmin(m) == 1){
        out(s,m) = log(cp * psimat(m,s) + DBL_MIN);
      } else if(Kmin(m) == 0){
        out(s,m) = log(cp * psimat(m,s) + (1-psimat(m,s)) + DBL_MIN);
      }
    }
  }

  return out;
}


// [[Rcpp::export]]
arma::mat get_loglik_pcount(arma::vec y, int M, arma::imat si,
                            arma::mat lammat, arma::mat pmat, int K,
                            arma::ivec Kmin){

  int S = lammat.n_cols;
  mat out(S,M);
  vec p_s(pmat.n_rows);
  vec psub, ysub;
  int J;
  double fac, ff, N, ky, numN;

  for (int s = 0; s < S; s++){

    p_s = pmat.col(s);

    for (int m = 0; m < M; m++){
      ysub = y.subvec(si(m,0), si(m,1));
      J = ysub.size();
      psub = p_s.subvec(si(m,0), si(m,1));
      fac = 1.0;
      ff = lammat(m,s) * prod(1 - psub);
      numN = K - Kmin(m);

      for (int i = 1; i < (numN + 1); i++){
        N = K - i + 1;
        ky = 1.0;
        for (int j = 0; j < J; j++){
          ky *= N / (N - ysub(j));
        }
        fac = 1 + fac * ff * ky / N;
      }

      out(s,m) = log(fac) + Rf_dpois(Kmin(m), lammat(m,s), true);
      for (int j = 0; j < J; j ++){
        out(s,m) += Rf_dbinom(ysub(j), Kmin(m), psub(j), true);
      }
    }
  }

  return out;
}
