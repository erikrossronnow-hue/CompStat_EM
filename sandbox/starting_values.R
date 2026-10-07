source("sandbox/setup.R")
library(ggplot2)
theme_set(theme_bw(base_size = 12))
set.seed(123)
nu = 2
y = 2 + 1.5 * rt(1000, df = nu)
starts = list("c(0, 1)" = c(0, 1), "c(6, 0.5)" = c(6, 0.5), "c(-4, 8)" = c(-4, 8))
out = list()
fits = list()

for (name in names(starts)) {
  tr = em_tracer(y, nu)
  fit = em_t_vec(starts[[name]], y, nu = nu, epsilon = 1e-16, cb = tr$tracer)
  trace = tr$trace()
  trace$start = name
  out[[name]] = trace
  fits[[name]] = data.frame(start = name, iterations = max(trace$n), mu = fit[1], sigma2 = fit[2], ll = loglik_t(fit, y, nu))
}
out = do.call(rbind, out)
fits = do.call(rbind, fits)
print(fits, row.names = FALSE)
write.csv(fits, "figures/starting_values.csv", row.names = FALSE)
paths = rbind(data.frame(out[c("n", "start")], par = "mu", val = out$mu), data.frame(out[c("n", "start")], par = "sigma2", val = out$sigma2))
ref = data.frame(par = c("mu", "sigma2"), val = c(fits$mu[1], fits$sigma2[1]))
p = ggplot(paths, aes(n, val, color = start)) + geom_hline(data = ref, aes(yintercept = val), linetype = 2, inherit.aes = FALSE) + geom_line() + facet_wrap(~ par, scales = "free_y") + labs(x = "EM iteration", y = "Parameter value", color = "Initial (mu, sigma2)", title = "Three starting points on the same sample; dashed lines mark the MLE")
ggsave("figures/starting_values.png", p, width = 9, height = 4.5, dpi = 160)
