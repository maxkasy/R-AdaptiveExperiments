library(furrr)

simulate_regret_reduction = function(A, B, k, assignment) {
  theta_draw = rbeta(rep(1, k), A, B)
  successes = rbinom(rep(1, k), as.integer(assignment), theta_draw)
  theta_bar_posterior = (A + successes) / (A + B + assignment)
  
  dstar = which.max(theta_bar_posterior)
  dstar_prior = which.max(A / (A + B))
  
  return(theta_draw[dstar] - theta_draw[dstar_prior])
}

expected_return = function (Y,
                            D,
                            #outcomes and treatments thus far
                            k,
                            #number of treatments
                            C = rep(0, k),
                            #vector of treatment cost
                            Nt,
                            RR = 50000) {
  A = 1 + tapply(Y, D, sum, default = 0)
  B = 1 + tapply(1 - Y, D, sum, default = 0)
  
  p = DtchoiceThompsonProbabilities(Y, D, k, C, RR)
  q = p * (1 - p)
  q = q / sum(q)
  
  n_floor = floor(Nt * q) #number of units assigned to each treatment, rounded down
  remainder = Nt * q - n_floor
  units_remaining = Nt - sum(n_floor) #remaining number of units to be assigned
  
  assignment = rep(0, k)
  if (units_remaining > 0)
    assignment[sample(1:k,
                      size = units_remaining,
                      replace = F,
                      prob = remainder)] = 1
  assignment = n_floor + assignment
  
  plan(multiprocess)
  regret_reduction_draws = future_map_dbl(1:RR, function(i)
    simulate_regret_reduction(A, B, k, assignment))
  
  return(mean(regret_reduction_draws))
}
