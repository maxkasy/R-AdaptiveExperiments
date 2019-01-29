betabinomial = function(n,s,a,b) {
  #probability mass function of the beta-binomial distribution
  #n trials, s successes, beta parameters a,b
  #in logs for numerical reasons
  exp(lgamma(n + 1)-lgamma(s + 1)-lgamma(n - s + 1) +
        lbeta((a + s),(b + n - s)) - lbeta(a,b)) 
}

# betabinomialvector = function(N,S,A,B) {
#   #wrapper for betabinomial:
#   #calculating betabinomial probability, product across components
#   prod(mapply(betabinomial, N,S,A,B))
# }


#requires global variable alpha0 specifying prior precision
betaposterior=function(D,Y){
  A=tapply(Y, D,sum) + alpha0 #posterior parameters, starting with Beta(alpha0,alpha0) (uniform) prior
  B=tapply(1-Y, D,sum) + alpha0
  list(A=A,B=B, theta= A /(A+B))
}

SWF= function(A, B, C){
  #social welfare function for expected average outcome of optimal treatment
  # A,B vectors with a, b parameters of beta distribution
  # C vector with cost of treatments
  max(A /(A+B) - C)
}


PolicyChoice=function(A,B,C){
  # policy choice maximizing expected average outcome of optimal treatment
  # A,B vectors with a, b parameters of beta distribution
  # C vector with cost of treatments
  which.max(A /(A+B) - C)
}



Regret=function(D, #vector of treatments across all waves
                Y, #vector of outcomes across all waves
                C, #treatment cost vector
                theta #true parameter vector
)
{
  k=max(D) #number of treatment arms, assuming we start at 1
  posterior=betaposterior(D,Y)
  d=PolicyChoice(posterior$A,posterior$B,C) # policy chosen given experimental data
  regret=max(theta-C)-theta[d]+C[d] # regret - difference in welfare between chosen and optimal option
  c(regret, d)
}


#dividing N units equally (up to rounding) between k treatments  
EqualAssignment=function(N,k){
  floor((1:N-1)*k/N)+1
}

#D vector based on vector n of assigned units
GivenAssignment=function(n,k){
  Dt=as.vector(unlist(sapply(1:k, function(i) rep(i,n[i])))) 
}



# beginning of period value function U
# expected welfare as a fuction of design n
# for prior A, B, , cost C, end of period value function Vfunction
U=function(A,B,C,n, Vfunction=SWF){
  k=length(A) #number of treatment arms
  
  #possible success vectors
  nplus1=n+1;
  N=prod(nplus1)
  SM=matrix(0, N, k)
  i=0:(N-1)
  for (j in k:1) { #mapping each element of i into corresponding vector of successes
    SM[,j]=i%%nplus1[j]
    i=i%/%nplus1[j]
  }
  
  #probabilities of successes for each treatment separately
  pp=map(1:k, function(j) betabinomial(n[j],0:n[j], A[j], B[j]))
  #probabilities of successes in SM matrix, stored as dataframe
  ppp=map_dfc(1:k, function(j) pp[[j]][SM[,j] + 1])
  #probability of each combination of successes, multiplying across treatment arms
  p= map_dbl(1:N, function(ii) prod(ppp[ii,]))
  
  #posterior expected social welfare for each combination of successes
  SW=sapply(1:N, function(ii) Vfunction(A+SM[ii,] ,B+n-SM[ii,] ,C))
  
  #expected SW
  p %*% SW
}


# for each design calculate U, given sample size N
# maximum over these will give value function V
UoverSimplex=function(A,B,C,N,
                      Ufunction=U, #how to estimate expected welfare (U or Uhat or Vfunction...)
                      coverage="full"){  #what assignments to try
  k=length(A)
   
  #number of units assigned to each of k treatments, summing to N
  # depending on value of coverage, do full, random, or handpicked
  nmatrix=simplex(N,k, coverage, RR=400, thetahat=A/(A+B)) 
  
  
  USimplex=as.data.frame(nmatrix)
  Nassignments=nrow(nmatrix) #number of different treatment assignments
  names(USimplex)=paste("n", 1:k, sep="")
  
  #calculate Ufunction for each row of the nmatrix
  USimplex$U=sapply(1:Nassignments, function(i) Ufunction(A,B,C, nmatrix[i,]))

  USimplex
}




# beginning of period value function V
# maximizing over possible designs n
# for prior A, B, , cost C
# sample sizes NN for coming waves
V=function(A,B,C,NN) {
  if (length(NN)>1) {
    Ufunction=function(A,B,C,n) U(A,B,C,n, Vfunction=function(A,B,C) V(A,B,C,NN[-1]))
  } else {
    Ufunction=U
  }
  
  USimplex=UoverSimplex(A,B,C,NN[1], Ufunction)
  max(USimplex$U) 
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
    k=length(A)
    thetadraws=sapply(1:k, function(j) rbeta(Nt, A[j], B[j]))
    Dt=sapply(1:Nt, function(j) which.max(thetadraws[j,]))
  } else if (method=="modifiedthompson") {
    # assign treatment in proportion to probability of being best
    # but don't allow same treatment two times in a row
    k=length(A)
    Dt=rep(0,Nt)
    Dt[1]=which.max(sapply(1:k, function(j) rbeta(1, A[j], B[j])))
    i=2
    while (i<=Nt) {
      Dt[i]=which.max(sapply(1:k, function(j) rbeta(1, A[j], B[j])))
      if (Dt[i] !=Dt[i-1]) {i=i+1}
    }
  }
  
  Dt
}
