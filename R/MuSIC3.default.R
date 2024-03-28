MuSIC3.default <- function(dataset_incident, dataset_prevalent, xw_rule, cl, optim.method = "nlminb", ...) {
  
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

  MLE_trunc <- optimf(c(-0.1, -0.1), likelihood_prevalent, par_coef = est0, ind_coef = ind_coef,
                  dataset_prevalent = dataset_prevalent, xw_rule = xw_rule, cl = NULL, csum = TRUE)
  est0_trunc <- MLE_trunc$par
  if (any(exp(est0_trunc) < 0.0001)) est0_trunc <- c(-0.1, -0.1)
  
  MLE_pall <- optimf(c(est0, est0_trunc), likelihood_forhess, likelihood_forhess_grad, ind_coef = ind_coef,
                     dataset_incident = dataset_incident, dataset_prevalent = dataset_prevalent, 
                     cl = NULL, csum = TRUE, xw_rule = xw_rule)
  est_all <- MLE_pall$par
  
  # Jackknife
  uID_inc <- unique(dataset_incident[,1])
  uID_pre <- unique(dataset_prevalent[,1])
  
  list_jack <- cbind(c(rep("inc", n_inc), rep("pre", n_pre)), c(uID_inc, uID_pre))
  
  fjack <- function(jack, list_jack, initial, dataset_inc, dataset_pre, xw_rule, ...) {
    
    uID_inc <- unique(dataset_inc[,1])
    uID_pre <- unique(dataset_pre[,1])
    
    ind_dat1 <- list_jack[jack, 1]
    ind_dat2 <- as.numeric(list_jack[jack, 2])
    if(ind_dat1 == "pre") subdat <- dataset_pre
    
    remained <- subdat[subdat[,1] != uID_inc[ind_dat2],]
    
    if(ind_dat1 == "inc") dataset_inc <- remained
    if(ind_dat1 == "pre") dataset_pre <- remained
    
    try_b <- try(res_jack <- optimf(initial, likelihood_forhess, likelihood_forhess_grad, ind_coef = ind_coef,
                                    dataset_incident = dataset_inc, dataset_prevalent = dataset_pre, 
                                    cl = NULL, csum = TRUE, xw_rule = xw_rule, control = list(trace = T)))
    
    res = list()
    if(class(try_b) != "try-error") {
      res$est_all <- res_jack$est_all
      res$est_inc <- res_jack$est_inc
    } else {
      res$est_all <- rep(NA, 21)
      res$est_inc <- rep(NA, 19)
    }
    
    res
  }
  
  njacks <- nrow(list_jack)
  res_jack <- mclapply(as.list(1:njacks), fjack, list_jack = list_jack, initial = est_all, dataset_inc = dataset_incident, dataset_pre = dataset_prevalent, 
                       xw_rule = xw_rule, mc.cores = length(cl)) 
  
  est_all_jacks <- matrix(NA, nrow = njacks, ncol = 21)
  est_inc_jacks <- matrix(NA, nrow = njacks, ncol = 19)
  
  for (jj in 1:njacks) {
    est_all_jacks[jj,] <- res_jack[[jj]]$est_all
    est_inc_jacks[jj,] <- res_jack[[jj]]$est_inc
  }
  
  results <- list()
  results$results_all <- results_all
  results$est_all_jacks <- est_all_jacks
  results$est_inc_jacks <- est_inc_jacks
  results$dataset_incident <- dataset_incident
  results$dataset_prevalent <- dataset_prevalent
  
  return(results)
  
}
