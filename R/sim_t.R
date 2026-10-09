sim_t = function(n, par, nu) {
  w = rgamma(n, shape = nu / 2, rate = 1 / 2)
  rnorm(n, mean = par[1], sd = sqrt(nu * par[2] / w))
}
