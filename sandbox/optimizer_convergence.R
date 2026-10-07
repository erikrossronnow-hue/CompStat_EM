source("sandbox/setup.R")
library(ggplot2)
theme_set(theme_bw(base_size = 12))
set.seed(123)
nu = 4
y = 2 + 1.5 * rt(50000, df = nu)
par = c(0, 1)
ref = em_t_cpp(par, y, nu = nu, epsilon = 1e-24)
ll_ref = loglik_t(ref, y, nu)
funs = list("EM naive" = em_t, "EM vectorized" = em_t_vec, "EM Rcpp steps" = em_t_cpp)
out = list()

# These traces include callback overhead; benchmarks.R times untraced fits.
for (name in names(funs)) {
  tr = em_tracer(y, nu)
  fit = funs[[name]](par, y, nu = nu, epsilon = 1e-20, cb = tr$tracer)
  trace = tr$trace()
  out[[name]] = data.frame(method = name, time = trace$time, ll = trace$ll)
}

for (method in c("BFGS", "nlminb")) {
  rows = list(data.frame(time = 0, ll = loglik_t(par, y, nu)))
  t0 = proc.time()[["elapsed"]]
  nll = function(q) {
    ll = loglik_t(c(q[1], exp(q[2])), y, nu)
    rows[[length(rows) + 1]] <<- data.frame(time = proc.time()[["elapsed"]] - t0, ll = ll)
    -ll / length(y)
  }
  grad = function(q) -colSums(score_t(c(q[1], exp(q[2])), y, nu)) * c(1, exp(q[2])) / length(y)
  q = c(par[1], log(par[2]))
  if (method == "BFGS") {
    fit = optim(q, nll, gr = grad, method = method, control = list(reltol = 1e-12, maxit = 10000))
  } else {
    fit = nlminb(q, nll, gradient = grad, control = list(rel.tol = 1e-12, sing.tol = 0, iter.max = 10000))
  }
  rows[[length(rows) + 1]] = data.frame(time = proc.time()[["elapsed"]] - t0, ll = loglik_t(c(fit$par[1], exp(fit$par[2])), y, nu))
  trace = do.call(rbind, rows)
  out[[method]] = data.frame(method = method, time = trace$time, ll = cummax(trace$ll))
}
out = do.call(rbind, out)
out$gap = pmax(ll_ref - out$ll, 1e-10)
write.csv(out, "figures/optimizer_convergence.csv", row.names = FALSE)
plotdat = rbind(transform(out, view = "All methods"), transform(subset(out, method != "EM naive"), view = "Faster methods"))
p = ggplot(plotdat, aes(time, gap, color = method)) + geom_step() + scale_y_log10() + facet_wrap(~ view, scales = "free_x") + labs(x = "Elapsed time (seconds)", y = "Best attained log-likelihood gap (floor at 1e-10)", color = NULL, title = "Convergence against real time on the same 50,000 observations", subtitle = "One traced run per method, including tracing overhead")
ggsave("figures/optimizer_convergence.png", p, width = 11, height = 5, dpi = 160)
