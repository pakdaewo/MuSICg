likelihood_forhess_DA <- function(est_all, ind_coef, dataset_incident, dataset_prevalent, xw_rule, cl, ...) {
  
  np1 <- ind_coef[length(ind_coef)] - 1
  
  par_coef <- est_all[1:np1]
  par_trunc <- est_all[-c(1:np1)]

  l_incident <- likelihood_incident_DA(par_coef, ind_coef, dataset_incident, cl = NULL)
  
  if (is.null(dataset_prevalent)) {
    l_prevalent <- 0
  } else {
    l_prevalent <- likelihood_prevalent_DA(par_trunc, par_coef, ind_coef, dataset_prevalent, xw_rule, cl) 
  }
  
  l_incident + l_prevalent
  
}

# # 
# l <- likelihood_forhess(c(par_coef, est0_trunc), ind_coef, dataset_incident, dataset_prevalent, xw_rule, cl)
# g <- likelihood_forhess_grad(c(par_coef, est0_trunc), ind_coef, dataset_incident, dataset_prevalent, xw_rule, cl)
# 
# gn <- numDeriv::grad(func = likelihood_forhess, x = c(par_coef, est0_trunc), ind_coef = ind_coef,
#                      dataset_incident = dataset_incident, dataset_prevalent = dataset_prevalent, xw_rule = xw_rule, cl = cl)
# 
# g - gn
# gn
