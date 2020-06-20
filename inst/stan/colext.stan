functions{

//can shortcut here I think
vector get_pY(int[] y, vector logit_p, int nd){
  vector[2] out;
  out[1] = nd;
  out[2] = exp(bernoulli_logit_lpmf(y | logit_p));
  return out;
}

matrix phi_matrix(vector phi_raw){
  return to_matrix(phi_raw, 2, 2, 0);
}

//delta-step transition prob matrix via Chapman-Kolmogorov equation
matrix get_phi(vector phi_raw, int Tstart, int Tnext){
  int delta = Tnext - Tstart;
  if(delta == 1){
    return phi_matrix(phi_raw[Tstart,]);
  }

  matrix phi = diag_matrix(rep_vector(1, 2));
  for (d in 1:delta){
    phi = phi * phi_matrix(phi_raw[(Tstart + d - 1),]);
  }
  return phi;
}

//Ts = indices of primary periods when site was sampled (eg not all NA)
real lp_colext(int[] y, int[] Ts, int[] J, vector psi, matrix phi_raw, 
               vector p, int[] nd){
  
  int T = size(Ts);
  matrix phi_prod = diag_matrix(rep_vector(1, 2));
  matrix phi;
  vector Dpt;
  int idx = 1;
  int end;

  for (t in 1:(T-1)){
    phi = get_phi(phi_raw, Ts[t], Ts[t+1]);
    end = idx + J[t] - 1;
    Dpt = get_pY(y[idx:end], p[idx:end], nd[t]);
    phi_prod *= diag_pre_multiply(Dpt, phi);
    idx += J[t];
  }

  end = idx + J[T] - 1;
  Dpt = get_pY(y[idx:end], p[idx:end], nd[T]);

  return log(dot_product(psi * phi_prod, Dpt));
}

//needs fixed
vector get_loglik_colext(int[] y, int M, int Tmax, int T, int[,] J, 
                         matrix psi_raw, matrix phi_raw, vector logit_p){
  vector[M] out;
  int idx = 1;
  int phi_idx = 1;
  int t_idx = 1;
  int end;
  int phi_end;
  int t_end;
  for (i in 1:M){
    end = idx + sum(J[,i]) - 1;
    phi_end = phi_idx + Tmax - 1;
    t_end = t_idx + T[i] - 1;

    out[i] = lp_colext(y[idx:end], Ts[t_idx:t_end], J[,i], psi[i,], 
                       phi_raw[phi_idx:phi_end,], logit_p[idx:end], nd[,i]);
    idx += sum(J[,i];
    phi_idx += Tmax;
    ts_idx += T[i];
  }
  return out;
}

}

data{

#include /include/data_single_season.stan

}

transformed data{

int no_detects[M];
for (m in 1:M){
  no_detects[m] = 1 - Kmin[m];
}

}

parameters{

#include /include/params_single_season.stan

}

transformed parameters{

vector[M] logit_psi;
vector[R] logit_p;
vector[M] log_lik;

logit_psi = X_state * beta_state;
logit_p = X_det * beta_det;

if(has_random_state){
  logit_psi = logit_psi + 
              csr_matrix_times_vector(Zdim_state[1], Zdim_state[2], Zw_state,
                                      Zv_state, Zu_state, b_state);
}
if(has_random_det){
  logit_p = logit_p + 
            csr_matrix_times_vector(Zdim_det[1], Zdim_det[2], Zw_det,
                                    Zv_det, Zu_det, b_det);
}

log_lik = get_loglik_occu(y, M, J, logit_psi, logit_p, no_detects);

}

model{

#include /include/model_single_season.stan

}

