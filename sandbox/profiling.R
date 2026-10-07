source("sandbox/setup.R")
library(ggplot2)
theme_set(theme_bw(base_size = 12))
dir.create("work/profiles", recursive = TRUE, showWarnings = FALSE)
set.seed(123)
nu = 4
y = 2 + 1.5 * rt(20000, df = nu)
par = c(0, 1)
funs = list(Naive = em_t, Vectorized = em_t_vec, Rcpp = em_t_cpp, Rcpp_ll = function(...) em_t_cpp(..., criterion = "ll"))
out = list()

# Repeat for two seconds so the profiler can sample these short calls.
for (name in names(funs)) {
  fit = funs[[name]]
  fit(par, y, nu = nu, epsilon = 1e-12)
  file = paste0("work/profiles/", name, ".out")
  p = profvis::profvis({
    t0 = proc.time()[["elapsed"]]
    while (proc.time()[["elapsed"]] - t0 < 2) {
      fit(par, y, nu = nu, epsilon = 1e-12)
    }
  }, interval = 0.01, prof_output = file)
  htmlwidgets::saveWidget(p, paste0("work/profiles/", name, ".html"), selfcontained = FALSE, libdir = "libs")
  s = summaryRprof(file)
  print(name)
  print(head(s$by.self, 8))
  vals = head(s$by.self, 8)
  out[[name]] = data.frame(version = name, fun = gsub('"', "", rownames(vals)), seconds = vals$self.time, percent = vals$self.pct)
}
out = do.call(rbind, out)
write.csv(out, "figures/profiling.csv", row.names = FALSE)
p = ggplot(out, aes(reorder(fun, percent), percent)) + geom_col(fill = "steelblue") + coord_flip() + facet_wrap(~ version, scales = "free_y") + labs(x = NULL, y = "Share of sampled self time (%)", title = "profvis: where each EM implementation spends its time", subtitle = "C++ appears as one call; its internal lines are not sampled")
ggsave("figures/profiling.png", p, width = 10, height = 7, dpi = 160)
