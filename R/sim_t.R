sim_t = function(n, par, v) {
  w = rgamma(n, shape = v / 2, rate = 1 / 2)
  rnorm(n, mean = par[1], sd = sqrt(v * par[2] / w))
}
