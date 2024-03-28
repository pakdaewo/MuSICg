likelihood_incident_DA <- function(par_coef, ind_coef, dataset_incident, cl, ...) {
  
  log_lambda <- par_coef[ind_coef[1]:(ind_coef[2] - 1)] # alpha_01, ..., alpha_23
  num_trans <- length(log_lambda)
  gamma_all <- exp(par_coef[ind_coef[2]:(ind_coef[3] - 1)]) # b_0, b_1, ..., b_nknots
  par_beta <- matrix(par_coef[ind_coef[3]:(ind_coef[4] - 1)], ncol = num_trans - 1, byrow = T) # x_1's ab and x_2's ab
  par_beta <- cbind(par_beta[,1], c(0, 0), par_beta[, -1])
  
  # Q (time-invariant variable X)
  ID <- dataset_incident[, 1]
  cen <- dataset_incident[, 2]
  visit <- dataset_incident[, 3]
  state <- dataset_incident[, 4]
  
  unique_ID <- unique(ID)
  len_ID <- length(unique_ID)
  X_ID <- dataset_incident[cumsum(table(ID)), 4 + 1:nrow(par_beta)] # bring X from the last observation at each ID
  
  Xb <- X_ID %*% par_beta
  alpha_ab <- matrix(log_lambda, nrow = 1)[rep(1, len_ID), ]
  Qe <- exp(alpha_ab + Xb)
  
  target_ID <- rbind(unique_ID, 1:length(unique_ID))
  l_target <- split(target_ID, rep(1:ncol(target_ID), each = nrow(target_ID)))
  
  #inherits(cl, "cluster")
  
  if (is.null(cl)) {
    log_trans_prob <- t(apply(target_ID, 2, fun_lik_incident_i, ID = ID, visit = visit, state = state, Qe = Qe, X_ID = X_ID, cen = cen, gamma_all = gamma_all))
  } else {
    if (.Platform$OS.type != "unix") {
      parallel_fun <- if (isTRUE(getOption("pboptions")$use_lb)) parLapplyLB else parLapply
      collected <- parallel_fun(cl = cl, l_target, fun_lik_incident_i, ID = ID, visit = visit, state = state, Qe = Qe, X_ID = X_ID, cen = cen, gamma_all = gamma_all)
    } else {
      collected <- mclapply(l_target, fun_lik_incident_i, ID = ID, visit = visit, state = state, Qe = Qe, X_ID = X_ID, cen = cen, gamma_all = gamma_all, 
                            mc.cores = length(cl))
    }
    log_trans_prob <- t(sapply(collected, function(x) x))
  }
  
  loglik <- sum(log_trans_prob) #- sum(par_knotsb^2)/(2 * par_theta^2)
  
  -loglik
}


# system.time(likelihood_incident(par_coef, ind_coef, dataset_incident, cl=NULL))
# system.time(likelihood_incident(par_coef, ind_coef, dataset_incident, cl=NULL))

# likelihood_incident(par_coef, ind_coef, dataset_incident, cl=NULL)
# # 
# g <- grad(func = likelihood_incident, x = par_coef, ind_coef = ind_coef, dataset_incident = dataset_incident, cl = NULL)
# gn <- likelihood_incident_grad(par_coef, ind_coef, dataset_incident, cl = NULL)
# g - gn

# system.time(gn <- likelihood_incident_grad(par_coef, ind_coef, dataset_incident, cl = NULL))
# system.time(gn <- likelihood_incident_grad(par_coef, ind_coef, dataset_incident, cl = cl))
