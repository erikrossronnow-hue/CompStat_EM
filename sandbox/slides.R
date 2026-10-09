source("sandbox/setup.R")
library(ggplot2)
theme_set(theme_bw(base_size = 16) + theme(legend.position = "bottom", panel.grid.minor = element_blank()))
b = read.csv("figures/benchmarks.csv")
b$criterion = factor(b$criterion, c("par", "ll"), c("Parameter change", "Likelihood change"))
p = ggplot(b, aes(n, ms, color = version)) + geom_line() + geom_point() + scale_x_log10() + scale_y_log10() + facet_wrap(~ criterion) + labs(x = "n", y = "Time [ms]", color = NULL)
ggsave("figures/benchmark_slides.png", p, width = 10, height = 4, dpi = 180)
pdat = subset(read.csv("figures/profiling.csv"), version == "Naive")
pdat$fun = factor(pdat$fun, c("<GC>", "m_step", "e_step"), c("GC", "M-step", "E-step"))
p = ggplot(pdat, aes(fun, percent)) + geom_col(fill = "steelblue") + coord_flip() + labs(x = NULL, y = "Self time [%]")
ggsave("figures/profile_slides.png", p, width = 4.5, height = 3.5, dpi = 180)
set.seed(123)
nu = 2
y = 2 + 1.5 * rt(5000, df = nu)
ref = em_t_vec(c(0, 1), y, nu = nu, epsilon = 1e-24)
rows = list()
for (criterion in c("par", "ll")) {
  for (eps in c(1e-6, 1e-10, 1e-14)) {
    tr = em_tracer(y, nu)
    fit = em_t_vec(c(0, 1), y, nu = nu, epsilon = eps, criterion = criterion, cb = tr$tracer)
    d = tr$trace()
    rows[[length(rows) + 1]] = data.frame(d, rule = criterion, eps = factor(eps), gap = pmax(loglik_t(ref, y, nu) - d$ll, 1e-12))
  }
}
d = do.call(rbind, rows)
d$rule = factor(d$rule, c("par", "ll"), c("Parameter change", "Likelihood change"))
ends = d[!duplicated(d[c("rule", "eps")], fromLast = TRUE), ]
p = ggplot(d, aes(n, gap, group = eps)) + geom_line(color = "grey60") + geom_point(data = ends, aes(color = eps), size = 3) + scale_y_log10() + facet_wrap(~ rule) + labs(x = "Iteration", y = "Likelihood gap", color = "epsilon")
ggsave("figures/convergence_slides.png", p, width = 10, height = 4, dpi = 180)
set.seed(123)
y = 2 + 1.5 * rt(1000, df = nu)
starts = list("(0, 1)" = c(0, 1), "(6, 0.5)" = c(6, 0.5), "(-4, 8)" = c(-4, 8))
rows = list()
for (name in names(starts)) {
  tr = em_tracer(y, nu)
  fit = em_t_vec(starts[[name]], y, nu = nu, epsilon = 1e-16, cb = tr$tracer)
  d = tr$trace()
  rows[[name]] = rbind(data.frame(n = d$n, start = name, par = "mu", val = d$mu), data.frame(n = d$n, start = name, par = "sigma2", val = d$sigma2))
}
d = do.call(rbind, rows)
p = ggplot(d, aes(n, val, color = start)) + geom_line() + facet_wrap(~ par, scales = "free_y") + labs(x = "Iteration", y = "Estimate", color = "Start")
ggsave("figures/starts_slides.png", p, width = 10, height = 4, dpi = 180)
d = subset(read.csv("figures/optimizer_convergence.csv"), method %in% c("EM Rcpp steps", "GD", "BFGS", "nlminb"))
d$method[d$method == "EM Rcpp steps"] = "EM"
p = ggplot(d, aes(time, gap, color = method)) + geom_step() + scale_y_log10() + labs(x = "Time [s]", y = "Likelihood gap", color = NULL)
ggsave("figures/optimizer_slides.png", p, width = 10, height = 4, dpi = 180)
set.seed(123)
nu = 4
y = sim_t(1000, c(2, 2.25), nu)
par = c(0, 1)
fits = rbind(EM = em_t_vec(par, y, nu = nu, epsilon = 1e-20), GD = fit_t_gd(par, y, nu)$par, BFGS = fit_t(par, y, nu)$par, nlminb = fit_t(par, y, nu, method = "nlminb")$par)
write.csv(data.frame(method = rownames(fits), fits, row.names = NULL), "figures/estimates_slides.csv", row.names = FALSE)
phat = fits[1, ]
info = list(Observed = obs_fisher_t(phat, y, nu), Empirical = emp_fisher_t(phat, y, nu), Expected = fisher_t(phat, length(y), nu), SEM = sem_fisher_t(phat, y, nu))
out = do.call(rbind, lapply(names(info), function(name) data.frame(method = name, i11 = info[[name]][1, 1], i12 = info[[name]][1, 2], i22 = info[[name]][2, 2], se_mu = sqrt(solve(info[[name]])[1, 1]), se_sigma2 = sqrt(solve(info[[name]])[2, 2]))))
write.csv(out, "figures/fisher_slides.csv", row.names = FALSE)
print(fits)
print(out, row.names = FALSE)
