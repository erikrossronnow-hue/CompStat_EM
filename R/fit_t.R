fit_t = function(par, y, v, method = "BFGS", epsilon = 1e-12) {
  q = c(par[1], log(par[2]))
  nll = function(q) neg_loglik_t(q, y, v)
  grad = function(q) grad_neg_loglik_t(q, y, v)
  if (method == "nlminb") {
    # Disable PORT's premature singular-convergence test.
    fit = nlminb(q, nll, gradient = grad, control = list(rel.tol = epsilon, sing.tol = 0, iter.max = 10000))
    val = fit$objective
  } else {
    fit = optim(q, nll, gr = grad, method = method, control = list(reltol = epsilon, maxit = 10000))
    val = fit$value
  }
  list(par = c(mu = unname(fit$par[1]), sigma2 = unname(exp(fit$par[2]))), value = val * length(y), convergence = fit$convergence)
}
