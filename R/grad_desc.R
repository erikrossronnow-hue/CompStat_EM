# Adapted from CompStat_stochgrad/R/stoch_grad.R.
grad_desc = function(par, N, gamma, epoch = batch, ..., maxit = 10000, epsilon = 1e-20, sampler = seq_len, cb = NULL) {
  if (is.function(gamma)) gamma = gamma(seq_len(maxit))
  gamma = rep_len(gamma, maxit)
  for (n in seq_len(maxit)) {
    par0 = par
    samp = sampler(N)
    par = epoch(par, samp, gamma[n], ...)
    if (!is.null(cb)) cb()
    done = sum((par - par0)^2) <= epsilon * (sum(par^2) + epsilon)
    if (done) return(list(par = par, iterations = n, convergence = 0))
  }
  list(par = par, iterations = maxit, convergence = 1)
}
