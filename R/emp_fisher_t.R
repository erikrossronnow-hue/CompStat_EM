emp_fisher_t = function(par, y, v, center = FALSE) {
  grad = score_t(par, y, v)
  if (center) grad = sweep(grad, 2, colMeans(grad))
  crossprod(grad)
}
