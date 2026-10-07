emp_fisher_t = function(par, y, nu, center = FALSE) {
  grad = score_t(par, y, nu)
  if (center) grad = sweep(grad, 2, colMeans(grad))
  crossprod(grad)
}
