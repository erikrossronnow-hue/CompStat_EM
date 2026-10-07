# Information for n observations, in coordinates (mu, sigma2).
fisher_t = function(par, n, nu) diag(c(n * (nu + 1) / ((nu + 3) * par[2]), n * nu / (2 * (nu + 3) * par[2]^2)))
