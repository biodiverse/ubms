//Basic data setup for single-season model
int model_code;
int M;
int T;
int Tsamp_size;
int Tsamp[Tsamp_size];
int R;
int J[M,T];
int y[R];
int si[M, 6];
int K;
int Kmin[M,T];
int y_dist;
int z_dist;
int n_aux1;
int n_aux2;
int n_aux3;
int aux1[n_aux1]; //Used for various auxiliary data
vector[n_aux2] aux2;
vector[n_aux3] aux3;

int has_random_state;
int has_random_det;
int n_obs_state;
int n_obs_det;
int n_fixed_state;
int n_fixed_det;
int n_group_vars_state;
int n_group_vars_det;
int n_random_state[has_random_state ? n_group_vars_state : 1];
int n_random_det[has_random_det ? n_group_vars_det: 1];
matrix[n_obs_state, n_fixed_state] X_state;
matrix[n_obs_det, n_fixed_det] X_det;
vector[n_obs_state] offset_state;
vector[n_obs_det] offset_det;

int Zdim_state[5];
vector[Zdim_state[3]] Zw_state;
int Zv_state[Zdim_state[4]];
int Zu_state[Zdim_state[5]];

int Zdim_det[5];
vector[Zdim_det[3]] Zw_det;
int Zv_det[Zdim_det[4]];
int Zu_det[Zdim_det[5]];

// Stuff for custom priors
int prior_dist_state[3];
int prior_dist_det[3];
int prior_dist_shape[3];
int prior_dist_scale[3];

matrix[3, (n_fixed_state+1)] prior_pars_state;
matrix[3, (n_fixed_det+1)] prior_pars_det;
matrix[3, 2] prior_pars_shape;
matrix[3, 2] prior_pars_scale;
