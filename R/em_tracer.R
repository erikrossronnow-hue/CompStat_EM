em_tracer = function(y, v) {
  rows = list()
  tracer = function() {
    env = parent.frame()
    par = env$par
    par0 = env$par0
    if (env$i == 1) rows[[1]] <<- data.frame(n = 0, time = 0, mu = par0[1], sigma2 = par0[2], ll = loglik_t(par0, y, v), par_norm_diff = NA)
    rows[[env$i + 1]] <<- data.frame(n = env$i, time = proc.time()[["elapsed"]] - env$t0, mu = par[1], sigma2 = par[2], ll = loglik_t(par, y, v), par_norm_diff = sqrt(sum((par - par0)^2)))
  }
  list(tracer = tracer, trace = function() do.call(rbind, rows))
}
