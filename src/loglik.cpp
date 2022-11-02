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

// [[Rcpp::export]]
arma::mat get_loglik_occuRN(arma::vec y, int M, arma::imat si,
                            arma::mat lammat, arma::mat rmat, int K,
                            arma::ivec Kmin){

  int S = lammat.n_cols;
  mat out(S,M);
  vec r_s(rmat.n_rows);
  vec ysub, qsub;
  int J;
  double f, g, p, lik;

  for (int s = 0; s < S; s++){

    r_s = rmat.col(s);

    for (int m = 0; m < M; m++){
      lik = 0.0;
      ysub = y.subvec(si(m,0), si(m,1));
      qsub = 1 - r_s.subvec(si(m,0), si(m,1));
      J = ysub.size();
      for (int k=Kmin(m); k < (K+1); k++){
        f = Rf_dpois(k, lammat(m, s), false);
        g = 0.0;
        for (int j = 0; j < J; j++){
          p = 1 - pow(qsub(j), k);
          g += Rf_dbinom(ysub(j), 1, p, true);
        }
        lik += f * exp(g);
      }
      out(s,m) = log(lik + DBL_MIN);
    }
  }

  return out;

}

// pi functions
arma::vec pi_removal(arma::vec p){
  int J = p.size();
  vec pi_out(J);
  pi_out(0) = p(0);
  for (int j = 1; j < J; j++){
    pi_out(j) = pi_out(j-1) / p(j-1) * (1-p(j-1)) * p(j);
  }
  return pi_out;
}

arma::vec pi_double(arma::vec p){
  vec pi_out(3);
  pi_out(0) = p(0) * (1 - p(1));
  pi_out(1) = p(1) * (1 - p(0));
  pi_out(2) = p(0) * p(1);
  return pi_out;
}

arma::vec pi_fun(int pi_type, arma::vec p, int J){
  vec out(J);
  if(pi_type == 0){
    out = pi_double(p);
  } else if(pi_type == 1){
    out = pi_removal(p);
  } else{
    Rcpp::stop("Invalid pi function");
  }
  return out;
}

// [[Rcpp::export]]
arma::mat get_loglik_multinomPois(arma::vec y, int M, arma::imat si,
                                  arma::mat lammat, arma::mat pmat, int pi_type){

  int S = lammat.n_cols;
  mat out = zeros(S,M);
  vec p_s(pmat.n_rows);
  vec ysub, psub;
  int R = pmat.n_rows / M;
  int pstart, pend, J;
  vec cp;

  for (int s = 0; s < S; s++){

    p_s = pmat.col(s);
    pstart = 0;

    for (int m = 0; m < M; m++){
      pend = pstart + R - 1;
      ysub = y.subvec(si(m,0), si(m,1));
      psub = p_s.subvec(pstart, pend);
      J = ysub.size();
      cp = pi_fun(pi_type, psub, J);

      for (int j = 0; j < J; j++){
        out(s, m) += Rf_dpois(ysub(j), lammat(m, s) * cp(j), true);
      }
      pstart += R;
    }
  }
  return out;
}


arma::vec get_pY(arma::vec y, arma::vec p, int nd){

  vec out(2);
  int J = y.size();
  out(0) = nd;
  out(1) = 1.0;
  for (int j = 0; j < J; j++){
    out(1) *= Rf_dbinom(y(j), 1, p(j), false);
  }
  return out;
}

arma::mat phi_matrix(arma::rowvec phi_raw){
  mat out(2,2);
  out(0,0) = phi_raw(0);
  out(0,1) = phi_raw(1);
  out(1,0) = phi_raw(2);
  out(1,1) = phi_raw(3);
  return out;
}

arma::mat get_phi(arma::mat phi_raw, int Tstart, int Tnext){
  mat phi = eye(2,2);
  int delta = Tnext - Tstart;
  if(delta == 1){
    return phi_matrix(phi_raw.row(Tstart));
  }

  for (int d = 1; d < (delta + 1); d++){
    phi = phi * phi_matrix(phi_raw.row(Tstart + d - 1));
  }
  return phi;
}



// [[Rcpp::export]]
arma::mat get_loglik_colext(arma::vec y, int M, arma::ivec Tsamp, arma::imat J,
                            arma::imat si, arma::cube psicube, arma::cube phicube,
                            arma::mat pmat, arma::imat nd){

  int S = pmat.n_cols;
  mat out(S,M);
  vec p_s(pmat.n_rows);
  mat psi_s(psicube.n_rows, psicube.n_cols);
  mat phi_s(phicube.n_rows, phicube.n_cols);

  rowvec psi_m;
  mat phi_m;
  ivec Tsamp_m;
  irowvec J_m, nd_m;
  vec y_m, p_m;
  int T;
  mat phi_prod(2,2);
  mat phi(2,2);
  vec Dpt(2);
  int idx, end;

  for (int s = 0; s < S; s++){

    p_s = pmat.col(s);
    psi_s = psicube.slice(s);
    phi_s = phicube.slice(s);

    for (int m = 0; m < M; m++){
      idx = 0;
      y_m = y.subvec(si(m,0), si(m,1));
      Tsamp_m = Tsamp.subvec(si(m,2), si(m,3));
      J_m = J.row(m);

      psi_m = psi_s.row(m);
      phi_m = phi_s.rows(si(m,4),si(m,5));
      p_m = p_s.subvec(si(m,0), si(m,1));
      nd_m = nd.row(m);

      T = Tsamp_m.size();
      phi_prod = eye(2,2);

      if(T > 1){
        for (int t = 0; t < (T-1); t++){
          phi = get_phi(phi_m, Tsamp_m(t), Tsamp_m(t+1));
          end = idx + J_m(Tsamp_m(t)) - 1;
          Dpt = get_pY(y_m.subvec(idx, end), p_m.subvec(idx, end),
                       nd_m(Tsamp_m(t)));

          phi_prod *= diagmat(Dpt) * phi;
          idx += J_m(Tsamp_m(t));
        }
      }

      end = idx + J_m(Tsamp_m(T-1)) - 1;
      Dpt = get_pY(y_m.subvec(idx,end), p_m.subvec(idx,end), nd_m(Tsamp_m(T-1)));
      out(s, m) = log(dot(psi_m * phi_prod, Dpt));
    }
  }

  return out;
}


arma::vec ttd_prob_exp(arma::vec y, arma::vec lam, arma::ivec delta){
  int J = y.size();
  vec e_lamt(J);

  for (int j=0; j < J; j++){
    e_lamt(j) = pow(lam(j), delta(j)) * exp(-lam(j)*y(j));
  }
  return e_lamt;
}

arma::vec ttd_prob_weib(arma::vec y, arma::vec lam, arma::ivec delta, double k){
  int J = y.size();
  vec e_lamt(J);
  for (int j = 0; j < J; j++){
    e_lamt(j) = pow(k*lam(j)*pow(lam(j)*y(j), (k-1)), delta(j)) *
                exp(-1*pow(lam(j)*y(j), k));
  }
  return e_lamt;
}


// [[Rcpp::export]]
arma::mat get_loglik_occuTTD(arma::vec y, int M, arma::imat si, arma::mat psimat,
                             arma::mat lammat, arma::vec k,
                             arma::ivec delta, int ydist){

  int S = psimat.n_cols;
  mat out(S,M);
  vec lam_s(lammat.n_rows);
  vec e_lamt;
  vec y_m, lam_m;
  ivec delta_m;
  double lik;

  for (int s = 0; s < S; s++){

    lam_s = lammat.col(s);

    for (int m = 0; m < M; m++){
      y_m = y.subvec(si(m,0), si(m,1));
      lam_m = lam_s.subvec(si(m,0), si(m,1));
      delta_m = delta.subvec(si(m,0), si(m,1));

      if(ydist == 1){
        e_lamt = ttd_prob_exp(y_m, lam_m, delta_m);
      } else if(ydist == 3){
        e_lamt = ttd_prob_weib(y_m, lam_m, delta_m, k(s));
      }

      lik = psimat(m,s) * prod(e_lamt) + (1-psimat(m,s)) * (1-max(delta_m));
      out(s,m) = log(lik);
    }
  }

  return out;
}
