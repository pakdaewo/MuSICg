likelihood_forhess_grad_DA <- function(est_all, ind_coef, dataset_incident, dataset_prevalent, xw_rule, cl, csum = TRUE, ...) {
  
  np1 <- ind_coef[length(ind_coef)] - 1
  
  par_coef <- est_all[1:np1]
  par_trunc <- est_all[-c(1:np1)]
  
  l_incident <- likelihood_incident_grad_DA(par_coef, ind_coef, dataset_incident, cl = cl, csum = csum)
  
  if (is.null(dataset_prevalent)) {
    l_prevalent <- 0
  } else {
    l_prevalent <- likelihood_prevalent_grad_DA(par_trunc, par_coef, ind_coef, dataset_prevalent, xw_rule, cl, csum = csum) 
  }
  
  if (csum == TRUE) lik =  c(l_incident,0,0) + l_prevalent
  if (csum == FALSE) {
    lik = rbind(cbind(l_incident, 0, 0), l_prevalent)
  } 
  
  return(lik)
  
}


# system.time(l_incident <- likelihood_incident_grad(par_coef, ind_coef, dataset_incident, cl = NULL, csum = TRUE))
# system.time(l_incident <- likelihood_incident_grad(par_coef, ind_coef, dataset_incident, cl = cl, csum = TRUE))
