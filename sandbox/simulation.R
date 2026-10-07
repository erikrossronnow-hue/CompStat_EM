source("sandbox/setup.R")
library(ggplot2)
theme_set(theme_bw(base_size = 12))
set.seed(123)
out = list()

for (nu in c(1, 4, 12)) {
  for (n in c(300, 1000)) {
    for (b in 1:40) {
      y = 2 + 1.5 * rt(n, df = nu)
      fit = em_t_vec(c(0, 1), y, nu = nu, epsilon = 1e-16)
      out[[length(out) + 1]] = data.frame(nu = nu, n = n, rep = b, mu = fit[1], sigma2 = fit[2])
    }
  }
}
out = do.call(rbind, out)
write.csv(out, "figures/simulation.csv", row.names = FALSE)
errs = rbind(data.frame(out[c("nu", "n", "rep")], par = "mu", error = out$mu - 2), data.frame(out[c("nu", "n", "rep")], par = "sigma2", error = out$sigma2 - 2.25))
print(aggregate(error ~ nu + n + par, errs, function(x) c(bias = mean(x), rmse = sqrt(mean(x^2)))))
p = ggplot(errs, aes(factor(nu), error, fill = factor(n))) + geom_hline(yintercept = 0, linetype = 2) + geom_boxplot(outlier.size = 1) + facet_wrap(~ par, scales = "free_y") + labs(x = "Degrees of freedom", y = "Estimate minus true parameter", fill = "Sample size", title = "EM estimates over 40 simulated samples per setting")
ggsave("figures/simulation.png", p, width = 9, height = 4.5, dpi = 160)

set.seed(123)
nu = 4
y = 2 + 1.5 * rt(1000, df = nu)
fit = em_t_vec(c(0, 1), y, nu = nu, epsilon = 1e-16)
x = seq(quantile(y, 0.01), quantile(y, 0.99), length.out = 400)
dens = rbind(data.frame(x = x, density = dt((x - 2) / 1.5, nu) / 1.5, curve = "True"), data.frame(x = x, density = dt((x - fit[1]) / sqrt(fit[2]), nu) / sqrt(fit[2]), curve = "Fitted"))
p = ggplot(data.frame(y = y), aes(y)) + geom_histogram(aes(y = after_stat(density)), bins = 60, fill = "grey85", color = "white") + geom_line(data = dens, aes(x, density, color = curve), linewidth = 0.8, inherit.aes = FALSE) + coord_cartesian(xlim = range(x)) + labs(x = "Observation", y = "Density", color = NULL, title = "Fixed nu = 4: fitted and true t density")
ggsave("figures/fitted_density.png", p, width = 8, height = 4.5, dpi = 160)
