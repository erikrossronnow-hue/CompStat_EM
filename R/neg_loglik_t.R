# q = c(mu, log(sigma2)), mean negative marginal log likelihood.
neg_loglik_t = function(q, y, v) -loglik_t(c(q[1], exp(q[2])), y, v) / length(y)
