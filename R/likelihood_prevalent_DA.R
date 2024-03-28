likelihood_prevalent_DA <- function(par_trunc, par_coef, ind_coef, dataset_prevalent, xw_rule, cl, ...) {

  eta <- exp(par_trunc[1])
  k <- exp(par_trunc[2])

  par_knotsa <- par_coef[ind_coef[1]:(ind_coef[2] - 1)] # alpha_01, ..., alpha_23
  num_trans <- length(par_knotsa)
  gamma_all <- exp(par_coef[ind_coef[2]:(ind_coef[3] - 1)]) # b_0, b_1, ..., b_nknots
  par_beta <- matrix(par_coef[ind_coef[3]:(ind_coef[4] - 1)], ncol = num_trans - 1, byrow = T) # x_1's ab and x_2's ab
  par_beta <- cbind(par_beta[,1], c(0, 0), par_beta[, -1])
  
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
    lognum_den <- apply(target, 1, fn_lik_prevalent_i, u_ID = u_ID, v_ID = v_ID, Qe = Qe, cen = cen, gamma_all = gamma_all, 
                 eta = eta, k = k, px = px, pw = pw, px21 = px2[,1], px22 = px2[,2], pw2 = pw2)
  } else {
    if (.Platform$OS.type != "unix") {
      parallel_fun <- if (isTRUE(getOption("pboptions")$use_lb)) parLapplyLB else parLapply
      collected <- parallel_fun(cl = cl, l_target, fn_lik_prevalent_i, u_ID = u_ID, v_ID = v_ID, Qe = Qe, cen = cen, gamma_all = gamma_all, 
                                eta = eta, k = k, px = px, pw = pw, px21 = px2[,1], px22 = px2[,2], pw2 = pw2)
    } else {
      collected <- mclapply(l_target,  fn_lik_prevalent_i, u_ID = u_ID, v_ID = v_ID, Qe = Qe, cen = cen, gamma_all = gamma_all, 
                            eta = eta, k = k, px = px, pw = pw, px21 = px2[,1], px22 = px2[,2], pw2 = pw2, mc.cores = length(cl))
    }
    lognum_den <- sapply(collected, function(x) x)
  }

  loglik_pre <- sum(lognum_den)

  -loglik_pre
  
}

# tt <- function(x, ...) {
#   par_coef <- x[1:17]
#   par_trunc <- x[18:19]
#   likelihood_prevalent_DA(par_trunc, par_coef, ind_coef, dataset_prevalent, xw_rule, cl)
# }
# # l <- likelihood_prevalent(par_trunc, par_coef, ind_coef, dataset_prevalent, xw_rule, cl)
# g <- likelihood_prevalent_grad_DA(par_trunc, par_coef, ind_coef, dataset_prevalent, xw_rule, cl)
# gn <- grad(func = likelihood_prevalent_DA, par_trunc, par_coef = par_coef, ind_coef = ind_coef, dataset_prevalent = dataset_prevalent,
#      xw_rule = xw_rule, cl = cl)

# g1 <- grad(tt, c(par_coef, par_trunc))
