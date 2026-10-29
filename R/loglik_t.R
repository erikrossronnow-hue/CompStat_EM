loglik_t = function(par, y, v) sum(dt((y - par[1]) / sqrt(par[2]), df = v, log = TRUE)) - length(y) * log(par[2]) / 2
