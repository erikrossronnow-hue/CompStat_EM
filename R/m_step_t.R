m_step_t = function(w, y, v) {
  sw = 0
  sy = 0
  for (i in seq_along(y)) {
    sw = sw + w[i]
    sy = sy + w[i] * y[i]
  }
  mu = sy / sw
  s = 0
  for (i in seq_along(y)) {
    s = s + w[i] * (y[i] - mu)^2
  }
  c(mu = mu, sigma2 = s / (length(y) * v))
}
