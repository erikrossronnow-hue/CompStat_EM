fit_t = function(par, y, nu, method = "BFGS", epsilon = 1e-12) {
  q = c(par[1], log(par[2]))
  # Mean loss keeps trial steps small.
  nll = function(q) -loglik_t(c(q[1], exp(q[2])), y, nu) / length(y)
  grad = function(q) -colSums(score_t(c(q[1], exp(q[2])), y, nu)) * c(1, exp(q[2])) / length(y)
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
