# par = c(mu, sigma2), w is the latent variable in the assignment.
e_step_t = function(par, y, nu) {
  w = numeric(length(y))
  for (i in seq_along(y)) {
    w[i] = (nu + 1) / (1 + (y[i] - par[1])^2 / (nu * par[2]))
  }
  w
}
