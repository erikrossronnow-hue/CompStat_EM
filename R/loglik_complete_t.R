# W ~ Gamma(nu / 2, rate = 1 / 2), X | W ~ N(mu, nu * sigma2 / W).
loglik_complete_t = function(par, y, w, nu) sum(dnorm(y, par[1], sqrt(nu * par[2] / w), log = TRUE) + dgamma(w, nu / 2, rate = 1 / 2, log = TRUE))
