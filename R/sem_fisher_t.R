sem_fisher_t = function(par, y, nu) {
  iX = numDeriv::hessian(neg_Q_t, par, pp = par, y = y, nu = nu)
  Phi = function(pp) m_step_t_vec(e_step_t_vec(pp, y, nu), y, nu)
  DPhi = numDeriv::jacobian(Phi, par)
  (diag(2) - t(DPhi)) %*% iX
}
