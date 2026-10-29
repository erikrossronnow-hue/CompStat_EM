sem_fisher_t = function(par, y, v) {
  iX = numDeriv::hessian(neg_Q_t, par, pp = par, y = y, v = v)
  Phi = function(pp) m_step_t_vec(e_step_t_vec(pp, y, v), y, v)
  DPhi = numDeriv::jacobian(Phi, par)
  (diag(2) - t(DPhi)) %*% iX
}
