source("sandbox/setup.R")
set.seed(42)

# Integrating out W recovers the t density and its conditional mean.
par = c(2, 2.25)
nu = 4
for (x in c(-2, 2, 5)) {
  joint = function(w) dnorm(x, par[1], sqrt(nu * par[2] / w)) * dgamma(w, nu / 2, rate = 1 / 2)
  marginal = dt((x - par[1]) / sqrt(par[2]), nu) / sqrt(par[2])
  dens = integrate(joint, 0, Inf)$value
  ew = integrate(function(w) w * joint(w), 0, Inf)$value / dens
  stopifnot(abs(dens - marginal) < 1e-8, abs(ew - e_step_t_vec(par, x, nu)) < 1e-7)
}

# The complete-data M-step agrees with direct numerical maximization.
w = rgamma(200, nu / 2, rate = 1 / 2)
y = rnorm(200, par[1], sqrt(nu * par[2] / w))
complete = fit_complete_t(y, w, nu)
nll = function(q) -loglik_complete_t(c(q[1], exp(q[2])), y, w, nu)
direct = optim(c(0, 0), nll, method = "BFGS", control = list(reltol = 1e-13))
stopifnot(max(abs(complete - c(direct$par[1], exp(direct$par[2])))) < 1e-5)

count = 0
for (nu in c(1, 4, 12)) {
  for (n in c(50, 300, 1000)) {
    for (seed in 1:3) {
      set.seed(seed)
      y = 2 + 1.5 * rt(n, nu)
      par = c(0, 1)
      w = e_step_t(par, y, nu)
      stopifnot(isTRUE(all.equal(w, e_step_t_vec(par, y, nu))), isTRUE(all.equal(w, e_step_t_cpp(par, y, nu))))
      stopifnot(isTRUE(all.equal(m_step_t(w, y, nu), m_step_t_vec(w, y, nu))), isTRUE(all.equal(m_step_t(w, y, nu), m_step_t_cpp(w, y, nu))))
      num_grad = numDeriv::grad(function(par) loglik_t(par, y, nu), par)
      stopifnot(max(abs(num_grad - colSums(score_t(par, y, nu)))) < 1e-5)
      for (criterion in c("par", "ll")) {
        eps = if (criterion == "par") 1e-20 else 1e-14
        traces = list()
        fits = list()
        for (i in 1:3) {
          tr = em_tracer(y, nu)
          fits[[i]] = list(em_t, em_t_vec, em_t_cpp)[[i]](par, y, nu = nu, epsilon = eps, criterion = criterion, cb = tr$tracer)
          traces[[i]] = tr$trace()
          stopifnot(all(diff(traces[[i]]$ll) >= -1e-8), nrow(traces[[i]]) < 10001)
        }
        stopifnot(max(abs(fits[[1]] - fits[[2]])) < 1e-7, max(abs(fits[[1]] - fits[[3]])) < 1e-7)
        stopifnot(length(unique(vapply(traces, nrow, integer(1)))) == 1)
        for (method in c("BFGS", "nlminb")) {
          direct = fit_t(par, y, nu, method = method)
          stopifnot(direct$convergence == 0, max(abs(fits[[2]] - direct$par)) < 1e-4, abs(loglik_t(fits[[2]], y, nu) + direct$value) < 1e-7)
        }
        count = count + 1
      }
    }
  }
}

# SEM gives observed information, rather than complete-data information.
phat = em_t_vec(c(0, 1), y, nu = nu, epsilon = 1e-20)
obs = obs_fisher_t(phat, y, nu)
sem = sem_fisher_t(phat, y, nu)
stopifnot(max(abs(obs - sem)) / max(abs(obs)) < 1e-6, max(abs(sem - t(sem))) < 1e-5)
stopifnot(all(eigen(obs)$values > 0), all(eigen(emp_fisher_t(phat, y, nu))$values > 0))

# Expected information agrees with integration of the individual scores.
par = c(2, 2.25)
for (nu in c(1, 4, 12)) {
  dens = function(x) dt((x - par[1]) / sqrt(par[2]), nu) / sqrt(par[2])
  i11 = integrate(function(x) score_t(par, x, nu)[, 1]^2 * dens(x), -Inf, Inf)$value
  i22 = integrate(function(x) score_t(par, x, nu)[, 2]^2 * dens(x), -Inf, Inf)$value
  i12 = integrate(function(x) score_t(par, x, nu)[, 1] * score_t(par, x, nu)[, 2] * dens(x), -Inf, Inf)$value
  stopifnot(max(abs(diag(fisher_t(par, 1, nu)) - c(i11, i22))) < 1e-7, abs(i12) < 1e-7)
}
cat("Passed density, complete-data, gradient,", count, "EM comparisons, monotonicity and Fisher-information checks.\n")
