e_step_t_vec = function(par, y, v) (v + 1) / (1 + (y - par[1])^2 / (v * par[2]))
