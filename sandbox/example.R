source("sandbox/setup.R")
set.seed(123)
v = 4
y = sim_t(1000, c(2, 2.25), v)
par = c(0, 1)

fits = rbind(naive = em_t(par, y, v = v, epsilon = 1e-16),
  vec = em_t_vec(par, y, v = v, epsilon = 1e-16),
  cpp = em_t_cpp(par, y, v = v, epsilon = 1e-16),
  GD = fit_t_gd(par, y, v)$par,
  BFGS = fit_t(par, y, v)$par,
  nlminb = fit_t(par, y, v, method = "nlminb")$par)
print(fits)

phat = em_t_vec(par, y, v = v, epsilon = 1e-20)

# sigma2 is the squared t scale, not Var(X).
info = list(observed = obs_fisher_t(phat, y, v), empirical = emp_fisher_t(phat, y, v), expected = fisher_t(phat, length(y), v), SEM = sem_fisher_t(phat, y, v))
print(info)
print(sqrt(diag(solve(info$observed))))
