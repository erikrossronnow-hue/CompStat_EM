fit_t_gd = function(par, y, v, gamma = 0.5, epsilon = 1e-20, maxit = 10000, cb = NULL) {
  q = c(par[1], log(par[2]))
  grad = function(q, i) grad_neg_loglik_t(q, y, v, i)
  fit = grad_desc(q, length(y), gamma, grad = grad, m = length(y), epsilon = epsilon, maxit = maxit, cb = cb)
  par = c(mu = unname(fit$par[1]), sigma2 = unname(exp(fit$par[2])))
  list(par = par, value = -loglik_t(par, y, v), iterations = fit$iterations, convergence = fit$convergence)
}
