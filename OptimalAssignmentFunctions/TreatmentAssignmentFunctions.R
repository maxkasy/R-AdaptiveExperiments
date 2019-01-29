#dividing N units equally (up to rounding) between k treatments  
EqualAssignment=function(N,k){
  floor((1:N-1)*k/N)+1
}

#D vector based on vector n of assigned units
GivenAssignment=function(n,k){
  sample(rep(1:k, n))
}

#D vector based on shares
ProportionalAssignment = function(Shares, Nt) {
  k=length(Shares)
  n_floor=floor(Nt*Shares) #number of units assigned to each treatment, rounded down
  remainder = Nt*Shares - n_floor
  units_remaining= Nt-sum(n_floor) #remaining number of units to be assigned

  if (units_remaining > 0) {
    Dt= c(rep(1:k, n_floor), # assigning each treatment acoording to rounded down average count
          sample(1:k,size=units_remaining,replace=T, prob=remainder)) #remaining units assigned randomly from remain replicate assignments
  } else {
    Dt= rep(1:k, n_floor)
  }
  sample(Dt)
}

Thompson=function(A,B, Nt) {
  k=length(A)
  thetadraws=sapply(1:k, function(j) rbeta(Nt, A[j], B[j]))
  sapply(1:Nt, function(j) which.max(thetadraws[j,]))
}

ThompsonProbabilities=function(A,B,RR=1000) {
  tabulate(Thompson(A,B,RR), nbins=length(A))/RR
}


# alternative ways to pick last-wave design
Dtchoice=function(A,B,C,Nt, method="optimal"){
  k=length(A)
  if (method=="optimal") {
    # find the assignment maximizing expected welfare
    USimplex=UoverSimplex(A,B,C,Nt, Ufunction=U, coverage="full")
    nt=USimplex[which.max(USimplex$U), 1:k]
    Dt=GivenAssignment(nt,k)
  } else if (method == "optimalhat") {
    # find the assignment maximizing simulated expected welfare
    Seed(A,B, Nt)
    USimplex=UoverSimplex(A,B,C,Nt, Ufunction=Uhat)
    nt=USimplex[which.max(USimplex$U), 1:k]
    Dt=GivenAssignment(nt,k)
  } else if (method == "optimalrandom"){
    USimplex=UoverSimplex(A,B,C,Nt, Ufunction=U, coverage="random")
    nt=USimplex[which.max(USimplex$U), 1:k]
    Dt=GivenAssignment(nt,k)
  } else if (method == "optimalhatrandom") {
    # find the assignment maximizing simulated expected welfare, over random set of assignments
    Seed(A,B, Nt)
    USimplex=UoverSimplex(A,B,C,Nt, Ufunction=Uhat, coverage="random")
    nt=USimplex[which.max(USimplex$U), 1:k]
    Dt=GivenAssignment(nt,k)    
  } else if (method == "optimalhandpicked"){
    USimplex=UoverSimplex(A,B,C,Nt, Ufunction=U, coverage="handpicked")
    nt=USimplex[which.max(USimplex$U), 1:k]
    Dt=GivenAssignment(nt,k)    
  } else if (method == "besthalf"){
    # find the assignment which implements the better half of treatments with equal shares
    k=length(A)
    thetahat=A/(A+B)
    #get ordered list of treatments, starting with the highest
    bestoptions=order(thetahat,decreasing=TRUE)
    k2=ceiling(k/2) #consider half the options in second round. maybe do something more sophisticated here?
    #start with nt vector as if options were ordered, then reorder afterwards
    nt=c(rep(floor(Nt/k2),k2), rep(0,k-k2))
    nt[seq_len(Nt-sum(nt))]=nt[1:(Nt-sum(nt))]+1
    nt[bestoptions]=nt
    Dt=GivenAssignment(nt,k)
  } else if (method=="thompson") {
    # assign treatment in proportion to probability of being best
    Dt=Thompson(A,B, Nt)
  } else if (method=="expectedthompson") {
      alpha=ThompsonProbabilities(A,B)
      Dt=ProportionalAssignment(alpha, Nt)    
  } else if (method=="modifiedthompson") {
    # assign treatment in proportion to probability of being best
    # but don't allow same treatment two times in a row
    k=length(A)
    RT=5 #number of replicate draws over which to average
    DtRT=rep(0,Nt*RT)
    
    previousD=-Inf # auxiliary variable to avoid repeat assignments of same D
    for (i  in 1:(Nt*RT)) {
      thetadraw=sapply(1:k, function(j) rbeta(1, A[j], B[j]))
      DtRT[i]=which.max(thetadraw)
      if (DtRT[i] == previousD) {
        thetadraw[previousD] = -Inf
        DtRT[i]=which.max(thetadraw)
      }
      previousD = DtRT[i]
    }
    
    Shares=tabulate(DtRT,k) / RT #average count for each treatment value and covariate value, replicated sample
    Dt=ProportionalAssignment(Shares, Nt)
  } else if (method=="toptwo") {
    beta=.5
    alpha=ThompsonProbabilities(A,B)
    if (max(alpha)<1) {
      ratio=alpha/(1-alpha)
      Shares=alpha*(beta + (1-beta)*(sum(ratio)-ratio))
    } else {
      Shares=alpha
    }
    Dt=ProportionalAssignment(Shares, Nt)
  } else if (method=="modifiedthompsonanalytic") {
    alpha=ThompsonProbabilities(A,B)
    if (max(alpha)<1) {
      Shares=alpha*(1-alpha) # Based on stationary distribution for modified Thompson
      Shares=Shares/sum(Shares)
    } else {
      Shares=alpha
    }  
    Dt=ProportionalAssignment(Shares, Nt)  
  } else if (method=="weightedmodifiedthompson") {
    alpha=ThompsonProbabilities(A,B)
    thetabar=A/(A+B)
    delta=max(thetabar)-thetabar
    dstar=which.max(thetabar)

    Shares=alpha*delta    
    if (sum(Shares)>0) {
      Shares=.5*Shares/sum(Shares)    
      Shares[dstar]=.5
    } else {
      Shares=alpha
    }      

    Dt=ProportionalAssignment(Shares, Nt)   
  }
  
  # else if (method=="marginalvalue") {
  #   Dt=MaxMarginalValue(A,B,C,Nt)
  # }
  
  Dt
}


# 
# MarginalValue=function(theta,C, nn){
#   value=theta-C
#   dstar=which.max(value)
#   Delta=value[dstar]-value
#   sigma=sqrt(theta*(1-theta)/nn)
#   sigma_plus_1=sqrt(theta*(1-theta)/(nn+1))
#   
#   ProbaMistake=pnorm(-Delta/sqrt(sigma^2 + sigma[dstar]^2))
#   ProbaMistake_plus_1=pnorm(-Delta/sqrt(sigma_plus_1^2 + sigma[dstar]^2))
#   
#   # This gets the marginal value for all sub-optimal treatments
#   MargVal=-Delta*(ProbaMistake_plus_1-ProbaMistake)
#   # Marginal value for the optimal treatment:
#   #MargVal[dstar]=-sum(Delta*(pnorm(-Delta/sqrt(sigma^2 + sigma_plus_1[dstar]^2))
#   #                          - ProbaMistake) )
#   MargVal
# }
# 
# MaxMarginalValue=function(A,B,C,Nt){
#   Dt=rep(0,Nt)
#   k=length(A)
#   nn=rep(0,k)
#   mm=A+B
#   
#   for (i in 1:Nt) {
#     theta=map2_dbl(A, B, function(a,b) rbeta(1,a,b))
#     if (sample(0:1,1)==1) {
#       MargVal=MarginalValue(theta,C, mm+nn)
#       Dt[i]=which.max(MargVal)
#     } else {
#       Dt[i]=which.max(theta-C)
#     }  
#     nn[Dt[i]]=nn[Dt[i]]+1
#   }
#   Dt
# }

