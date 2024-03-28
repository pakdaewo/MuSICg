MuSIC2.default <- function(dataset_incident, dataset_prevalent, xw_rule, cl = NULL, optim.method = "nlminb", hess = FALSE, ...) {
  
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
  
  MLE_p <- optimf(par_coef, likelihood_incident, gradient = likelihood_incident_grad, ind_coef = ind_coef,
                  dataset_incident = dataset_incident, cl = NULL, csum = TRUE, control = list(maxit = 10000, trace = T))
  
  est0 <- est_inc <- MLE_p$par
  
  if (hess == TRUE) {
    hessian_inc <- numDeriv::hessian(func = likelihood_incident, x = est0, ind_coef = ind_coef, 
                                     dataset_incident = dataset_incident, cl = NULL)
  } else {
    hessian_inc <- NULL
  }
  
  len_par0 <- length(est0)
  crossprod_inc_i <- likelihood_incident_grad(est0, ind_coef, dataset_incident, cl, csum = FALSE)
  crossg_inc <- matrix(fn_sumxtx(crossprod_inc_i), len_par0, len_par0)
  
  # ## step 2 - finding truncation parameters
  # MLE_trunc <- optimf(c(0,0), likelihood_prevalent, gradient = likelihood_prevalent_grad, 
  #                     par_coef = est0, ind_coef = ind_coef, dataset_prevalent = dataset_prevalent, 
  #                     xw_rule = xw_rule, cl = cl, only.trunc = TRUE, control = list(maxit = 10000, trace = T))
  # 
  # est0_trunc <- MLE_trunc$par
  # 
  MLE_pall <- optimf(c(est0, -0.1, -0.1), likelihood_forhess, gradient = likelihood_forhess_grad, ind_coef = ind_coef,
                   dataset_incident = dataset_incident, dataset_prevalent = dataset_prevalent, 
                   cl = cl, csum = TRUE, xw_rule = xw_rule, control = list(maxit = 10000, trace = T))
  est_all <- MLE_pall$par
  
  # hessian
  if (hess == TRUE) {
    hessian_h <- numDeriv::hessian(func = likelihood_forhess, x = est_all, ind_coef = ind_coef,
                                   dataset_incident = dataset_incident, dataset_prevalent = dataset_prevalent,
                                   dist_trunc = dist_trunc, xw_rule = xw_rule, cl = cl)
  } else {
    hessian_h <- NULL
  }
  
  # crossproduct of gradient
  len_par <- length(est_all)
  cross_i <- likelihood_forhess_grad(est_all, ind_coef, dataset_incident, dataset_prevalent, xw_rule, cl, csum = FALSE)
  crossprod_g <- matrix(fn_sumxtx(cross_i), len_par, len_par)
  
  results <- list()
  results$est_all <- est_all
  results$hessian_h <- hessian_h
  results$crossprod_g <- crossprod_g
  results$est_inc <- est_inc
  results$hessian_inc <- hessian_inc
  results$crossg_inc <- crossg_inc
  results$dataset_incident <- dataset_incident
  results$dataset_prevalent <- dataset_prevalent
  
  return(results)
  
}
