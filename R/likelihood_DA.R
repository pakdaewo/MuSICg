likelihood_DA <- function(par_coef, ind_coef, par_trunc, dataset_incident, dataset_prevalent, xw_rule, cl, ...) {

  l_incident <- likelihood_incident_DA(par_coef, ind_coef, dataset_incident, cl = NULL)

  if (is.null(dataset_prevalent)) {
    l_prevalent <- 0
  } else {
    l_prevalent <- likelihood_prevalent_DA(par_trunc, par_coef, ind_coef, dataset_prevalent, xw_rule, cl) 
  }

  l_incident + l_prevalent

}

# #
# l <- likelihood_DA(par_coef, ind_coef, est0_trunc, dataset_incident, dataset_prevalent, xw_rule, cl)
# g <- likelihood_grad_DA(par_coef, ind_coef, est0_trunc, dataset_incident, dataset_prevalent, xw_rule, cl)
# gn <- numDeriv::grad(func = likelihood_DA, x = par_coef, ind_coef = ind_coef, par_trunc = est0_trunc,
#                      dataset_incident = dataset_incident, dataset_prevalent = dataset_prevalent, xw_rule = xw_rule, cl = cl)
# 
# g - gn
# gn