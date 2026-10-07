# Terms depending on par only; pp stays fixed during differentiation.
neg_Q_t = function(par, pp, y, nu) {
  w = e_step_t_vec(pp, y, nu)
  length(y) * log(par[2]) / 2 + sum(w * (y - par[1])^2) / (2 * nu * par[2])
}
