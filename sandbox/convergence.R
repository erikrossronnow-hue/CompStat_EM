source("sandbox/setup.R")
library(ggplot2)
theme_set(theme_bw(base_size = 12))
set.seed(123)
nu = 2
y = 2 + 1.5 * rt(5000, df = nu)
par = c(0, 1)
ref = em_t_vec(par, y, nu = nu, epsilon = 1e-24)
ll_ref = loglik_t(ref, y, nu)
out = list()
fits = list()

# Equal epsilon has different accuracy under the two stopping rules.
for (criterion in c("par", "ll")) {
  for (eps in c(1e-6, 1e-10, 1e-14)) {
    tr = em_tracer(y, nu)
    fit = em_t_vec(par, y, nu = nu, epsilon = eps, criterion = criterion, cb = tr$tracer)
    trace = tr$trace()
    trace$criterion = criterion
    trace$epsilon = factor(eps)
    trace$gap = pmax(ll_ref - trace$ll, 1e-12)
    out[[length(out) + 1]] = trace
    fits[[length(fits) + 1]] = data.frame(criterion = criterion, epsilon = eps, iterations = max(trace$n), mu = fit[1], sigma2 = fit[2], ll_gap = ll_ref - loglik_t(fit, y, nu), par_error = sqrt(sum((fit - ref)^2)))
  }
}
out = do.call(rbind, out)
fits = do.call(rbind, fits)
print(fits, row.names = FALSE)
write.csv(fits, "figures/convergence.csv", row.names = FALSE)
ends = out[!duplicated(out[c("criterion", "epsilon")], fromLast = TRUE), ]
p = ggplot(out, aes(n, gap, color = epsilon, group = epsilon)) + geom_line(color = "grey60") + geom_point(data = ends, size = 3) + scale_y_log10() + facet_wrap(~ criterion, labeller = as_labeller(c(par = "Parameter change (book)", ll = "Likelihood change"))) + labs(x = "EM iteration", y = "Log-likelihood gap (floor at 1e-12)", color = "Epsilon", title = "Two stopping rules at three tolerances", subtitle = "Colored points show where each tolerance stops")
ggsave("figures/convergence.png", p, width = 9, height = 4.5, dpi = 160)

tr = em_tracer(y, nu)
fit = em_t_vec(par, y, nu = nu, epsilon = 1e-24, cb = tr$tracer)
trace = subset(tr$trace(), n > 0 & par_norm_diff > 0)
tail_trace = subset(trace, n > max(n) / 2)
rate = exp(coef(lm(log(par_norm_diff) ~ n, data = tail_trace))["n"])
Phi = function(pp) m_step_t_vec(e_step_t_vec(pp, y, nu), y, nu)
radius = max(abs(eigen(numDeriv::jacobian(Phi, fit))$values))
print(c(trace_rate = rate, jacobian_rate = radius))
p = ggplot(trace, aes(n, par_norm_diff)) + geom_line() + geom_point(size = 1) + scale_y_log10() + labs(x = "EM iteration", y = "Norm of parameter change", title = paste("Local EM convergence rate:", round(rate, 3)))
ggsave("figures/convergence_rate.png", p, width = 7, height = 4.5, dpi = 160)
