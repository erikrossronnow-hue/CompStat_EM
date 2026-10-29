source("sandbox/setup.R")
set.seed(123)
nu = 4
y = sim_t(1000, c(2, 2.25), nu)
par = c(0, 1)

fits = rbind(naive = em_t(par, y, nu = nu, epsilon = 1e-16),
  vec = em_t_vec(par, y, nu = nu, epsilon = 1e-16),
  cpp = em_t_cpp(par, y, nu = nu, epsilon = 1e-16),
  GD = fit_t_gd(par, y, nu)$par,
  BFGS = fit_t(par, y, nu)$par,
  nlminb = fit_t(par, y, nu, method = "nlminb")$par)
print(fits)

phat = em_t_vec(par, y, nu = nu, epsilon = 1e-20)

# sigma2 is the squared t scale, not Var(X).
info = list(observed = obs_fisher_t(phat, y, nu), empirical = emp_fisher_t(phat, y, nu), expected = fisher_t(phat, length(y), nu), SEM = sem_fisher_t(phat, y, nu))
print(info)
print(sqrt(diag(solve(info$observed))))
