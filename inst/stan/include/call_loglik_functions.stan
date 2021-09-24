if(model_code == 0){
  log_lik = get_loglik_occu(y, M, J, si, lp_state, lp_det, Kmin[,1]);
} else if(model_code == 1){
  log_lik = get_loglik_rn(y, M, J, si, lp_state, lp_det, K, Kmin[,1]);
} else if(model_code == 2){
  log_lik = get_loglik_pcount(y, M, J, si, lp_state, lp_det, z_dist,
                              log_scale, K, Kmin[,1]);
} else if(model_code == 4){
  log_lik = get_loglik_distsamp(y, M, aux2, si, lp_state, lp_det, z_dist,
                                log_scale, aux1[1], y_dist, aux3);
} else if(model_code == 5){
  log_lik = get_loglik_multinomPois(y, M, si, lp_state, lp_det, y_dist);
} else if(model_code == 6){
  log_lik = get_loglik_occuTTD(aux2, M, si, lp_state, lp_det, log_shape,
                               aux1, y_dist);
}
