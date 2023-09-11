//Basic data setup for single-season model
int model_code;
int M;
int T;
int Tsamp_size;
array[Tsamp_size] int Tsamp;
int R;
array[M,T] int J;
array[R] int y;
array[M, 6] int si;
int K;
array[M,T] int Kmin;
int y_dist;
int z_dist;
int n_aux1;
int n_aux2;
int n_aux3;
array[n_aux1] int aux1; //Used for various auxiliary data
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
array[has_random_state ? n_group_vars_state : 1] int n_random_state;
array[has_random_det ? n_group_vars_det: 1] int n_random_det;
matrix[n_obs_state, n_fixed_state] X_state;
matrix[n_obs_det, n_fixed_det] X_det;
vector[n_obs_state] offset_state;
vector[n_obs_det] offset_det;

array[5] int Zdim_state;
vector[Zdim_state[3]] Zw_state;
array[Zdim_state[4]] int Zv_state;
array[Zdim_state[5]] int Zu_state;

array[5] int Zdim_det;
vector[Zdim_det[3]] Zw_det;
array[Zdim_det[4]] int Zv_det;
array[Zdim_det[5]] int Zu_det;

// Stuff for custom priors
array[3] int prior_dist_state;
array[3] int prior_dist_det;
array[3] int prior_dist_shape;
array[3] int prior_dist_scale;

matrix[3, (n_fixed_state+1)] prior_pars_state;
matrix[3, (n_fixed_det+1)] prior_pars_det;
matrix[3, 2] prior_pars_shape;
matrix[3, 2] prior_pars_scale;
