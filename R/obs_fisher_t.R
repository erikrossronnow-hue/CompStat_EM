obs_fisher_t = function(par, y, nu) numDeriv::hessian(function(par) -loglik_t(par, y, nu), par)
