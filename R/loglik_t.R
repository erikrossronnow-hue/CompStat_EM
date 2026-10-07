loglik_t = function(par, y, nu) sum(dt((y - par[1]) / sqrt(par[2]), df = nu, log = TRUE)) - length(y) * log(par[2]) / 2
