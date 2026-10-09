grad_neg_loglik_t = function(q, y, nu, i = seq_along(y)) {
  s = exp(q[2])
  -colMeans(score_t(c(q[1], s), y[i], nu)) * c(1, s)
}
