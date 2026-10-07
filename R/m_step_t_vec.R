m_step_t_vec = function(w, y, nu) {
  mu = sum(w * y) / sum(w)
  c(mu = mu, sigma2 = sum(w * (y - mu)^2) / (length(y) * nu))
}
