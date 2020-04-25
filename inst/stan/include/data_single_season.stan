//Basic data setup for single-season model
int M;
int J[M];
int y[sum(J)];
int has_random_state;
int has_random_det;
int n_fixed_state;
int n_fixed_det;
int n_group_vars_state;
int n_group_vars_det;
int n_random_state[has_random_state ? n_group_vars_state : 1];
int n_random_det[has_random_det ? n_group_vars_det: 1];
matrix[M, n_fixed_state] X_state;
matrix[sum(J), n_fixed_det] X_det;

int Zdim_state[5];
vector[Zdim_state[3]] Zw_state;
int Zv_state[Zdim_state[4]];
int Zu_state[Zdim_state[5]];

int Zdim_det[5];
vector[Zdim_det[3]] Zw_det;
int Zv_det[Zdim_det[4]];
int Zu_det[Zdim_det[5]];
