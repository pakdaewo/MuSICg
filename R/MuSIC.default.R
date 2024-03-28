MuSIC.default <- function(dataset_incident, dataset_prevalent, xw_rule, cl = NULL, optim.method = "nlminb", ...) {

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
  
  hessian_inc <- numDeriv::hessian(func = likelihood_incident, x = est0, ind_coef = ind_coef, 
                                   dataset_incident = dataset_incident, cl = NULL)
  len_par0 <- length(est0)
  crossprod_inc_i <- likelihood_incident_grad(est0, ind_coef, dataset_incident, cl, csum = FALSE)
  crossg_inc <- matrix(fn_sumxtx(crossprod_inc_i), len_par0, len_par0)
  
  ## step 2 - finding truncation parameters
  MLE_trunc <- optimf(c(0,0), likelihood_prevalent, gradient = likelihood_prevalent_grad, 
                      par_coef = est0, ind_coef = ind_coef, dataset_prevalent = dataset_prevalent, 
                      xw_rule = xw_rule, cl = cl, only.trunc = TRUE, control = list(maxit = 10000, trace = T))
  
  est0_trunc <- MLE_trunc$par

  ## step 3 - iterations
  # given truncation, optimize betas from likelihood
  distance <- 1
  tol <- 0.01

  old_est0 <- est0
  old_trunc <- est0_trunc
  old_all <- c(old_est0, old_trunc)
  track <- old_all
  ntry <- 1
  
  while ((distance > tol) & (ntry < 100)) {

    MLE_all <- optimf(old_est0, likelihood, gradient = likelihood_grad, ind_coef = ind_coef,
                      par_trunc = old_trunc, dataset_incident = dataset_incident, dataset_prevalent = dataset_prevalent, 
                      xw_rule = xw_rule, cl = cl, control = list(maxit = 10000, trace = T))
    
    new_est0 <- MLE_all$par

    MLE_trunc <- optimf(old_trunc, likelihood_prevalent, gradient = likelihood_prevalent_grad, 
                        par_coef = new_est0, ind_coef = ind_coef, dataset_prevalent = dataset_prevalent, 
                        xw_rule = xw_rule, cl = cl, only.trunc = TRUE, control = list(maxit = 10000, trace = T))
    
    
    new_trunc <- MLE_trunc$par

    new_all <- c(new_est0, new_trunc)
    distance <- sqrt(crossprod(new_all - old_all)/crossprod(old_all))

    track <- rbind(track, new_all)
    old_est0 <- new_est0
    old_trunc <- new_trunc
    old_all <- c(old_est0, old_trunc)
    
    ntry <- ntry + 1
    print(distance)
    print(track)
    
  }

  # hessian
  hessian_h <- numDeriv::hessian(func = likelihood_forhess, x = old_all, ind_coef = ind_coef,
                                 dataset_incident = dataset_incident, dataset_prevalent = dataset_prevalent,
                                 dist_trunc = dist_trunc, xw_rule = xw_rule, cl = cl)

  # crossproduct of gradient
  len_par <- length(old_all)
  cross_i <- likelihood_forhess_grad(old_all, ind_coef, dataset_incident, dataset_prevalent, xw_rule, cl, csum = FALSE)
  crossprod_g <- matrix(fn_sumxtx(cross_i), len_par, len_par)
  
  results <- list()
  results$est_all <- old_all
  results$est <- old_est0
  results$est_trunc <- old_trunc
  results$track <- track
  results$ntry <- ntry
  results$hessian_h <- hessian_h
  results$crossprod_g <- crossprod_g
  results$distance <- distance
  results$est_inc <- est_inc
  results$hessian_inc <- hessian_inc
  results$crossg_inc <- crossg_inc
  results$dataset_incident <- dataset_incident
  results$dataset_prevalent <- dataset_prevalent

  return(results)

}
