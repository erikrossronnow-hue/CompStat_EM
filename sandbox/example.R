source("sandbox/setup.R")
set.seed(123)
nu = 4
y = 2 + 1.5 * rt(1000, df = nu)
par = c(0, 1)

fits = rbind(naive = em_t(par, y, nu = nu, epsilon = 1e-16), vec = em_t_vec(par, y, nu = nu, epsilon = 1e-16), cpp = em_t_cpp(par, y, nu = nu, epsilon = 1e-16), BFGS = fit_t(par, y, nu)$par, nlminb = fit_t(par, y, nu, method = "nlminb")$par)
print(fits)

tr = em_tracer(y, nu)
phat = em_t_vec(par, y, nu = nu, epsilon = 1e-20, cb = tr$tracer)
em_trace = tr$trace()
tail_trace = subset(em_trace, n > max(n) / 2 & par_norm_diff > 0)
log_rate_fit = lm(log(par_norm_diff) ~ n, data = tail_trace)
print(exp(coef(log_rate_fit)["n"]))

# sigma2 is the squared t scale, not Var(X).
info = list(observed = obs_fisher_t(phat, y, nu), empirical = emp_fisher_t(phat, y, nu), expected = fisher_t(phat, length(y), nu), SEM = sem_fisher_t(phat, y, nu))
print(info)
print(sqrt(diag(solve(info$observed))))
