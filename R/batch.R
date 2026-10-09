# Adapted from CompStat_stochgrad/R/batch.R.
batch = function(par, samp, gamma, grad, m = 50, ...) {
  M = floor(length(samp) / m)
  for (j in seq_len(M)) {
    i = samp[((j - 1) * m + 1):(j * m)]
    par = par - gamma * grad(par, i, ...)
  }
  par
}
