likelihood_grad <- function(par_coef, ind_coef, par_trunc, dataset_incident, dataset_prevalent, xw_rule, cl, csum = TRUE, ...) {
  
  len_par <- length(par_coef)
  l_incident <- likelihood_incident_grad(par_coef, ind_coef, dataset_incident, cl = cl, csum = csum)
  
  if (is.null(dataset_prevalent)) {
    l_prevalent <- 0
  } else {
    l_prevalent <- likelihood_prevalent_grad(par_trunc, par_coef, ind_coef, dataset_prevalent, xw_rule, cl, csum = csum) 
  }
  
  if (csum == TRUE) lik = l_incident + l_prevalent[1:len_par]
  if (csum == FALSE) lik = rbind(l_incident, l_prevalent[,1:len_par])
  
  return(lik)
  
}

