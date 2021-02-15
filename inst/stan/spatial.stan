functions{

#include /include/functions_occu.stan

// Source for this function:
// Clark A, Altwegg R. 2019. Efficient Bayesian analysis of occupancy models
// with logit link functions. Ecology and Evolution 9: 756–768.
real theta_lpdf(vector theta, real tau, matrix Qalpha, int n_eigen){
  return 0.5*(n_eigen*log(tau) - tau*quad_form(Qalpha, theta));
}

//real lp_occu_probit(int[] y, real raw_psi, vector raw_p, int Kmin){
//  real out;
//  real psi = Phi(raw_psi);
//  int J = num_elements(raw_p);
//  vector[J] p;
//  for (j in 1:J){
//    p[j] = Phi(raw_p[j]);
//  }
//  out = exp(bernoulli_lpmf(y | p)) * psi + (1 - Kmin) * (1-psi);
//  return log(out);
//}

//vector get_loglik_occu_probit(int[] y, int M, int[,] J, int[,] si, vector raw_psi,
//                  vector raw_p, int[] Kmin){
//  vector[M] out;
//  for (i in 1:M){
//    out[i] = lp_occu_probit(y[si[i,1]:si[i,2]], raw_psi[i],
//                     raw_p[si[i,1]:si[i,2]], Kmin[i]);
//   }
//  return out;
//}
//
}

data{
#include /include/data.stan
int n_aug_sites;
int n_eigen;
matrix[M+n_aug_sites,n_eigen] Kmat;
matrix[n_eigen,n_eigen] Qalpha;
matrix[n_aug_sites,n_fixed_state] X_aug;
}

transformed data{
//vector[n_eigen] zeros;
matrix[M+n_aug_sites,n_fixed_state] X_state_all;
X_state_all = append_row(X_state, X_aug);
//zeros = rep_vector(0, n_eigen);

}

parameters{
#include /include/params_single_season.stan
real<lower=0> tau;
}

transformed parameters{

vector[M+n_aug_sites] lp_state;
vector[n_obs_det] lp_det;
vector[M] log_lik;

lp_state = X_state_all * beta_state;
lp_det = X_det * beta_det;

lp_state += Kmat * b_state;

if(has_random_det){
  lp_det = lp_det +
            csr_matrix_times_vector(Zdim_det[1], Zdim_det[2], Zw_det,
                                    Zv_det, Zu_det, b_det);
}

if(model_code == 0){
  log_lik = get_loglik_occu(y, M, J, si, lp_state, lp_det, Kmin[,1]);
  //log_lik = get_loglik_occu_probit(y, M, J, si, lp_state, lp_det, Kmin[,1]);
}

}

model{

#include /include/rand_priors_single_season.stan
#include /include/fixed_priors_single_season.stan

tau ~ gamma(0.5, 0.005);
b_state ~ theta_lpdf(tau, Qalpha, n_eigen);
//b_state ~ multi_normal_prec(zeros, tau * Qalpha);

target += sum(log_lik);

}
