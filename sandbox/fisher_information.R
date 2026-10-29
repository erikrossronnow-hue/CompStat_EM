source("sandbox/setup.R")
set.seed(123)
v = 4
y = 2 + 1.5 * rt(1000, df = v)
phat = em_t_vec(c(0, 1), y, v = v, epsilon = 1e-20)
info = list(observed = obs_fisher_t(phat, y, v), empirical = emp_fisher_t(phat, y, v), expected = fisher_t(phat, length(y), v), SEM = sem_fisher_t(phat, y, v))
out = do.call(rbind, lapply(names(info), function(name) data.frame(method = name, i11 = info[[name]][1, 1], i12 = info[[name]][1, 2], i22 = info[[name]][2, 2], se_mu = sqrt(solve(info[[name]])[1, 1]), se_sigma2 = sqrt(solve(info[[name]])[2, 2]))))
print(out, row.names = FALSE)
write.csv(out, "figures/fisher_information.csv", row.names = FALSE)
