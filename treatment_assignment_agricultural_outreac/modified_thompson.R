#D vector based on shares
ProportionalAssignment = function(Shares, Nt) {
  k = length(Shares)
  n_floor = floor(Nt * Shares) #number of units assigned to each treatment, rounded down
  remainder = Nt * Shares - n_floor
  units_remaining = Nt - sum(n_floor) #remaining number of units to be assigned
  
  if (units_remaining > 0) {
    Dt = c(rep(1:k, n_floor),
           # assigning each treatment acoording to rounded down average count
           sample(
             1:k,
             size = units_remaining,
             replace = T,
             prob = remainder
           )) #remaining units assigned randomly from remain replicate assignments
  } else {
    Dt = rep(1:k, n_floor)
  }
  sample(Dt)
}


DtchoiceThompson = function(Y, D, #outcomes and treatments thus far
                            k, #number of treatments
                            C, #vector of treatment cost
                            Nt) {
  # number of observations for period t
  
  A = 1 + tapply(Y, D, sum, default = 0) 
  B = 1 + tapply(1-Y, D, sum, default = 0) 
  
  Dt = rep(0, Nt)
  previousD = -Inf # auxiliary variable to avoid repeat assignments of same D
  
  for (i  in 1:Nt) {
    thetadraw = sapply(1:k, function(j)
      rbeta(1, A[j], B[j]))
    Dt[i] = which.max(thetadraw - C)
    previousD = Dt[i]
  }
  
  factor(Dt, levels = 1:k)
}



DtchoiceThompsonProbabilities = function(Y, D, #outcomes and treatments thus far
                                         k, #number of treatments
                                         C = rep(0, k), #vector of treatment cost
                                         RR = 50000) {
  #number of replication draws
  
  # Repeat Thompson sampling RR times
  DtRR = DtchoiceThompson(Y, D, k, C, RR)
  P_Dt = table(DtRR) / RR #average count for each treatment value and covariate value, replicated sample
  
  P_Dt = as.tibble(matrix(P_Dt, 1, k))
  colnames(P_Dt) = paste(1:k)
  P_Dt
}


DtchoiceThompson_modified = function(Y, D, #outcomes and treatments thus far
                                     k, #number of treatments
                                     C = rep(0, k), #vector of treatment cost
                                     Nt,
                                     return_probabilities=F) {
  # number of observations for period t
  alpha = DtchoiceThompsonProbabilities(Y, D, k, C)
  if (max(alpha) < 1) {
    Shares = alpha * (1 - alpha) # Based on stationary distribution for modified Thompson
    Shares = Shares / sum(Shares)
  } else {
    Shares = alpha
  }
  
  if (return_probabilities) return(Shares)
    else return(ProportionalAssignment(Shares, Nt))
}
