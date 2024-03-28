MuSIC4.default <- function(dataset_incident, dataset_prevalent, xw_rule, cl, optim.method = "nlminb", hess = FALSE, g.dist, tau = NULL, ...) {
  
  match.call()
  
  if (!is.null(cl)) {
    if (.Platform$OS.type == "windows") {
      if (!inherits(cl, "cluster"))
        cl <- NULL
    } else {
      if (inherits(cl, "cluster")) {
        if (length(cl) < 2L)
          
          cl <- NULL
      } else {
        if (cl < 2)
          cl <- NULL
      }
    }
  }
  
  visit <- dataset_incident[, 3]
  nX <- ncol(dataset_incident) - 4
  
  ## step 1 - initial beta's
  num_trans <- 6
  par_knotsa <- rep(0, num_trans) # a_(01), ..., a_(23)
  par_knotsa[1] <- par_knotsa[1] - 0.01
  par_knotsb <- 0
  par_beta <- matrix(0, ncol = num_trans, nrow = nX)
  num_beta <- length(par_beta)
  
  # for data analysis
  # par_beta <- par_beta[, -2]
  # num_beta <- length(par_beta)
  # 
  par_coef <- c(par_knotsa, par_knotsb, par_beta)
  ind_coef <- cumsum(c(1, num_trans, 1, num_beta))
  
  if (optim.method == "nlminb") optimf <- nlminb
  if (optim.method == "optim") optimf <- optim
  
  MLE_p <- optimf(par_coef, likelihood_incident, likelihood_incident_grad, ind_coef = ind_coef,
                  dataset_incident = dataset_incident, cl = NULL, csum = TRUE)
  est0 <- est_inc <- MLE_p$par

  if (g.dist == "weibull") {
    MLE_trunc <- optimf(c(-0.1, -0.1), likelihood_prevalent, par_coef = est0, ind_coef = ind_coef,
                        dataset_prevalent = dataset_prevalent, xw_rule = xw_rule, cl = NULL, csum = TRUE)
    est0_trunc <- MLE_trunc$par
    if (any(exp(est0_trunc) < 0.0001)) est0_trunc <- c(-0.1, -0.1)
  }

  if (g.dist == "exp") {
    MLE_trunc <- optimf(0, likelihood_prevalent_exp, par_coef = est0, ind_coef = ind_coef,
                        dataset_prevalent = dataset_prevalent, xw_rule = xw_rule, cl = NULL, csum = TRUE)
    est0_trunc <- MLE_trunc$par
    if (any(exp(est0_trunc) < 0.0001)) est0_trunc <- c(-0.1)
  }

  if (g.dist == "length") {
    est0_trunc <- NULL
    if(is.null(tau)) tau <- max(dataset_prevalent[,3] + dataset_prevalent[,5])
  }
  
  MLE_pall <- optimf(c(est0, est0_trunc), likelihood_forhess, ind_coef = ind_coef,
                     dataset_incident = dataset_incident, dataset_prevalent = dataset_prevalent, 
                     cl = NULL, csum = TRUE, xw_rule = xw_rule, g.dist = g.dist, tau = tau)
  est_all <- MLE_pall$par

  # hessian
  if (hess == TRUE) {
    hessian_h <- numDeriv::hessian(func = likelihood_forhess, x = est_all, ind_coef = ind_coef,
                                   dataset_incident = dataset_incident, dataset_prevalent = dataset_prevalent,
                                   dist_trunc = dist_trunc, xw_rule = xw_rule, cl = cl, g.dist = g.dist, tau = tau)
  } else {
    hessian_h <- NULL
  }
  
  results <- list()
  results$est_all <- est_all
  results$hessian_h <- hessian_h

  return(results)
  
}
