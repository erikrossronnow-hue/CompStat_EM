source("sandbox/setup.R")
library(ggplot2)
theme_set(theme_bw(base_size = 12))
out = list()
versions = c("Naive", "Vectorized", "Rcpp steps")
nu = 4
par = c(0, 1)
eps = 1e-12

# Same data, start, epsilon and stopping rule, compilation is untimed.
for (n in c(500, 5000, 50000)) {
  set.seed(123)
  y = 2 + 1.5 * rt(n, df = nu)
  for (criterion in c("par", "ll")) {
    b = bench::mark(naive = em_t(par, y, nu = nu, epsilon = eps, criterion = criterion), vec = em_t_vec(par, y, nu = nu, epsilon = eps, criterion = criterion), cpp = em_t_cpp(par, y, nu = nu, epsilon = eps, criterion = criterion), iterations = 7, check = TRUE, filter_gc = FALSE)
    out[[length(out) + 1]] = data.frame(n = n, criterion = criterion, version = versions, ms = as.numeric(b$median) * 1000, mb = as.numeric(b$mem_alloc) / 1024^2)
    print(out[[length(out)]], row.names = FALSE)
  }
}
out = do.call(rbind, out)
write.csv(out, "figures/benchmarks.csv", row.names = FALSE)
p = ggplot(out, aes(n, ms, color = version)) + geom_line() + geom_point() + scale_x_log10() + scale_y_log10() + facet_wrap(~ criterion, labeller = as_labeller(c(par = "Parameter change (book)", ll = "Likelihood change"))) + labs(x = "Sample size", y = "Median complete-fit time (ms)", color = NULL, title = "bench::mark: identical convergence rule within each panel", subtitle = "7 measured fits per version; epsilon = 1e-12; compilation excluded")
ggsave("figures/benchmarks.png", p, width = 9, height = 4.5, dpi = 160)

# Optimizers have their own stopping rules; report achieved accuracy too.
set.seed(123)
y = 2 + 1.5 * rt(5000, df = nu)
ref = em_t_cpp(par, y, nu = nu, epsilon = 1e-24)
fits = list(EM = em_t_cpp(par, y, nu = nu, epsilon = 1e-20), GD = fit_t_gd(par, y, nu)$par, BFGS = fit_t(par, y, nu)$par, nlminb = fit_t(par, y, nu, method = "nlminb")$par)
b = bench::mark(EM = em_t_cpp(par, y, nu = nu, epsilon = 1e-20), GD = fit_t_gd(par, y, nu), BFGS = fit_t(par, y, nu), nlminb = fit_t(par, y, nu, method = "nlminb"), iterations = 15, check = FALSE, filter_gc = FALSE)
methods = data.frame(method = names(fits), ms = as.numeric(b$median) * 1000, ll_gap = vapply(fits, function(fit) loglik_t(ref, y, nu) - loglik_t(fit, y, nu), numeric(1)), par_error = vapply(fits, function(fit) sqrt(sum((fit - ref)^2)), numeric(1)))
print(methods, row.names = FALSE)
write.csv(methods, "figures/optimizer_benchmarks.csv", row.names = FALSE)
p = ggplot(methods, aes(method, ms)) + geom_col(fill = "steelblue") + labs(x = NULL, y = "Median complete-fit time (ms)", title = "EM versus marginal-likelihood optimization", subtitle = "Different internal stopping rules; attained likelihood gaps in optimizer_benchmarks.csv")
ggsave("figures/optimizer_benchmarks.png", p, width = 8, height = 4.5, dpi = 160)
saveRDS(sessionInfo(), "work/benchmark_session.rds")
