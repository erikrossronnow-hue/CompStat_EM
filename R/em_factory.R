# Adapted from https://cswr.nrhstat.org/em, with two stopping rules.
em_factory = function(e_step, m_step, eps = 1e-6, ll = NULL) {
  force(e_step); force(m_step); force(eps); force(ll)
  function(par, y, epsilon = eps, criterion = "par", cb = NULL, maxit = 10000, ...) {
    if (!is.null(cb)) t0 = proc.time()[["elapsed"]]
    if (criterion == "ll") ll0 = ll(par, y, ...)
    for (i in seq_len(maxit)) {
      par0 = par
      par = m_step(e_step(par, y, ...), y, ...)
      if (criterion == "par") {
        done = sum((par - par0)^2) <= epsilon * (sum(par^2) + epsilon)
      } else {
        ll1 = ll(par, y, ...)
        done = abs(ll1 - ll0) <= epsilon * (abs(ll1) + epsilon)
        ll0 = ll1
      }
      if (!is.null(cb)) cb()
      if (done) return(par)
    }
    stop("EM reached maxit")
  }
}
