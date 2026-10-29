# Information for n observations, in coordinates (mu, sigma2).
fisher_t = function(par, n, v) diag(c(n * (v + 1) / ((v + 3) * par[2]), n * v / (2 * (v + 3) * par[2]^2)))
