likelihood_prevalent_length <- function(par_coef, ind_coef, dataset_prevalent, xw_rule, cl, tau, ...) {
  
  log_lambda <- par_coef[ind_coef[1]:(ind_coef[2] - 1)] # alpha_01, ..., alpha_23
  num_trans <- length(log_lambda)
  gamma_all <- exp(par_coef[ind_coef[2]:(ind_coef[3] - 1)]) # b_0, b_1, ..., b_nknots
  par_beta <- matrix(par_coef[ind_coef[3]:(ind_coef[4] - 1)], ncol = num_trans, byrow = T) # x_1's ab and x_2's ab
  
  # Q (time-invariant variable X)
  cen <- dataset_prevalent[, 2]
  u_ID <- dataset_prevalent[, 3]
  v_ID <- dataset_prevalent[, 5]
  X_ID <- dataset_prevalent[, 6:7]
  n_ID <- length(u_ID)
  
  Xb <- X_ID %*% par_beta
  
  alpha_ab <- matrix(log_lambda, nrow = 1)[rep(1, n_ID), ]
  Qe <- exp(alpha_ab + Xb)
  
  px <- xw_rule[, 1]
  pw <- xw_rule[, 2]
  px2 <- as.matrix(expand.grid(px, px))
  pw2 <- apply(expand.grid(pw, pw), 1, prod)
  
  target <- matrix(1:n_ID, ncol = 1)
  l_target <- as.list(target)
  
  if (is.null(cl)) {
    # lognum_den <- apply(target, 1, fn_lik_prevalent_i, u_ID = u_ID, v_ID = v_ID, Qe = Qe, cen = cen, gamma_all = gamma_all, 
    #                     eta = eta, k = k, px = px, pw = pw, px21 = px2[,1], px22 = px2[,2], pw2 = pw2)
    loglik_pre <- fn_lik_prevalent_length(u_ID, v_ID, Qe, cen, gamma_all, tau, px, pw, px2[,1], px2[,2], pw2)
  } else {
    if (.Platform$OS.type != "unix") {
      parallel_fun <- if (isTRUE(getOption("pboptions")$use_lb)) parLapplyLB else parLapply
      collected <- parallel_fun(cl = cl, l_target, fn_lik_prevalent_i_length, u_ID = u_ID, v_ID = v_ID, Qe = Qe, cen = cen, gamma_all = gamma_all, 
                                tau = tau, px = px, pw = pw, px21 = px2[,1], px22 = px2[,2], pw2 = pw2)
    } else {
      collected <- mclapply(l_target,  fn_lik_prevalent_i_length, u_ID = u_ID, v_ID = v_ID, Qe = Qe, cen = cen, gamma_all = gamma_all, 
                            tau = tau, px = px, pw = pw, px21 = px2[,1], px22 = px2[,2], pw2 = pw2, mc.cores = length(cl))
    }
    lognum_den <- sapply(collected, function(x) x)
    loglik_pre <- sum(lognum_den)
  }
  
  
  -loglik_pre
  
}
