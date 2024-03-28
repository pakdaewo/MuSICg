likelihood_prevalent_grad <- function(par_trunc, par_coef, ind_coef, dataset_prevalent, xw_rule, cl, csum = TRUE, only.trunc = FALSE, ...) {
  
  eta <- exp(par_trunc[1])
  k <- exp(par_trunc[2])

  par_knotsa <- par_coef[ind_coef[1]:(ind_coef[2] - 1)] # alpha_01, ..., alpha_23
  num_trans <- length(par_knotsa)
  gamma_all <- exp(par_coef[ind_coef[2]:(ind_coef[3] - 1)]) # b_0, b_1, ..., b_nknots
  par_beta <- matrix(par_coef[ind_coef[3]:(ind_coef[4] - 1)], ncol = num_trans, byrow = T) # x_1's ab and x_2's ab
  
  len_par <- length(c(par_knotsa, gamma_all, par_beta, par_trunc))
  
  # Q (time-invariant variable X)
  cen <- dataset_prevalent[, 2]
  u_ID <- dataset_prevalent[, 3]
  v_ID <- dataset_prevalent[, 5]
  X_ID <- dataset_prevalent[, 6:7]
  n_ID <- length(u_ID)
  
  Xb <- X_ID %*% par_beta
  
  alpha_ab <- matrix(par_knotsa, nrow = 1)[rep(1, n_ID), ]
  Qe <- exp(alpha_ab + Xb)
  
  lam23 <- Qe[,6]
  
  px <- xw_rule[, 1]
  pw <- xw_rule[, 2]
  px2 <- as.matrix(expand.grid(px, px))
  pw2 <- apply(expand.grid(pw, pw), 1, prod)
  
  target <- matrix(1:n_ID, ncol = 1)
  l_target <- as.list(target)

  if (is.null(cl)) {
    # log_trans_prob <- t(apply(target, 1, fn_lik_prevalent_grad_i, u_ID = u_ID, v_ID = v_ID, X_ID = X_ID, Qe = Qe, cen = cen, gamma_all = gamma_all, 
    #                         eta = eta, k = k, px = px, pw = pw, px21 = px2[,1], px22 = px2[,2], pw2 = pw2, len_par = len_par))
    log_trans_prob <- fn_lik_prevalent_grad(u_ID, v_ID, X_ID, Qe, cen, gamma_all, eta, k, px, pw, px2[,1], px2[,2], pw2, len_par)
  } else {
    if (.Platform$OS.type != "unix") {
      parallel_fun <- if (isTRUE(getOption("pboptions")$use_lb)) parLapplyLB else parLapply
      collected <- parallel_fun(cl = cl, l_target, fn_lik_prevalent_grad_i, u_ID = u_ID, v_ID = v_ID, X_ID = X_ID, Qe = Qe, cen = cen, gamma_all = gamma_all, 
                                eta = eta, k = k, px = px, pw = pw, px21 = px2[,1], px22 = px2[,2], pw2 = pw2, len_par = len_par)
    } else {
      collected <- mclapply(l_target, fn_lik_prevalent_grad_i, u_ID = u_ID, v_ID = v_ID, X_ID = X_ID, Qe = Qe, cen = cen, gamma_all = gamma_all, 
                            eta = eta, k = k, px = px, pw = pw, px21 = px2[,1], px22 = px2[,2], pw2 = pw2, len_par = len_par, mc.cores = length(cl))
    }
    log_trans_prob <- t(sapply(collected, function(x) x))
  }
  
  log_trans_prob[,7] <- log_trans_prob[,7] * gamma_all
  log_trans_prob[,len_par-1] <- log_trans_prob[,len_par-1] * eta
  log_trans_prob[,len_par] <- log_trans_prob[,len_par] * k
  
  if (only.trunc == FALSE) {
    if (csum == TRUE)  d_loglik <- colSums(log_trans_prob)
    if (csum == FALSE) d_loglik <- log_trans_prob
  }
  
  if (only.trunc == TRUE) {
    if (csum == TRUE)  d_loglik <- colSums(log_trans_prob)[(len_par - 1):len_par]
    if (csum == FALSE) d_loglik <- log_trans_prob[, (len_par - 1):len_par]
  }
  return(-d_loglik)
  
}
  
# system.time(log_trans_prob <- t(apply(target, 1, fn_lik_prevalent_grad_i, u_ID = u_ID, v_ID = v_ID, X_ID = X_ID, Qe = Qe, cen = cen, gamma_all = gamma_all, 
#                           eta = eta, k = k, px = px, pw = pw, px21 = px2[,1], px22 = px2[,2], pw2 = pw2, len_par = len_par)))
# 
# 
# system.time(b <- fn_lik_prevalent_grad(u_ID, v_ID, X_ID, Qe, cen, gamma_all, eta, k, px, pw, px2[,1], px2[,2], pw2, len_par))
# colSums(log_trans_prob) - colSums(b)
