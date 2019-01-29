library(tidyverse)


# number of replication draws for Uhat
RU=1000


# function for drawing outcomes from sampling distribution
simulatedSample=function(D, #treatment vector
                         theta #success probabilities
)
{rbinom(length(D), 1, theta[D]) #bernoulli draws with probability theta(D)
}



# generating random draws once and for all
# using global assignment operator <<- 
# for estimation of Uhat
Seed = function(A,B, Nmax){
  set.seed(12231983)
  k=length(A)
  #draws from prior distribution of theta - dimensions (RU, k)
  theta <<- sapply(1:k, function(d)  rbeta(RU, A[d], B[d]))
  #draws of potential outcomes - dimensions (Nmax, RU, k)
  Yd <<- sapply(theta, function(y) rbernoulli(Nmax, p = y))
  Yd <<- array(Yd, c(Nmax, RU,k))
}


# a simulated approximation of U
Uhat=function(A,B,C,n, Vfunction=SWF){
  k=length(A)
  # number of successes for design n, for each of r=1..RU replicates
  SR=sapply(1:k, function(d) if (n[d]>1) {colSums(Yd[1:n[d],,d],1)}
                                else if (n[d]==1) {Yd[1,,d]}
                                else {rep(0,RU)} )
  # posterior expected social welfare for each realization of SR
  SW=sapply(1:RU, function(r) Vfunction(A+SR[r,] ,B+n-SR[r,] ,C))
  mean(SW)
}




#creating nmatrix that has all elements of simplex
simplex = function(N,k, 
                   coverage="full", #what part of the simplex to trace
                   RR=500, #draws, if coverage = "random"
                   thetahat=NULL # if coverage = "handpicked"
                   ){
  Nplus1=N+1
  Nplus1tok=Nplus1^k
  if (coverage=="full"){ 
    #computationally costly version of getting full simplex
    nmatrix=matrix(0, Nplus1tok, k)
    i=1:(Nplus1tok)
    for (j in k:1) {
      nmatrix[,j]=i%%Nplus1
      i=i%/%Nplus1
    }
    #keep only rows with correct rowsum
    nmatrix[rowSums(nmatrix)==N,,drop=FALSE]
  } else if (coverage=="random") { 
    #getting random subsample of simplex
    nmatrix=matrix(runif(RR*k), RR, k) #sample from contiuous hypercube
    nmatrix=t(apply(nmatrix, 1, function(x) N*x/sum(x))) #row normalize. this gives more points at center!
    nmatrix=floor(nmatrix) #round down
    nmatrix+  #add missing numbers to get right rowsum again
      t(sapply(N-rowSums(floor(nmatrix)), function(dN) sample(c(rep(1, dN), rep(0, k-dN))))) 
  } else if (coverage=="handpicked") { 
    # candidate set for optimal assignments 
    #get ordered list of treatments, starting with the highest
    bestoptions=order(thetahat,decreasing=TRUE)
    # candidates: equal treatment between j best options, and between rest
    
    # auxiliary function to divide up sample
     divider=function(j,nn){ #j best treatments preferred, nn units assigned to each of those
      ndiv=c(rep(nn,j), rep(floor((N-j*nn)/(k-j)),k-j))
      remainder=N-sum(ndiv)
      ndiv[seq_len(remainder)]=ndiv[seq_len(remainder)]+1
      ndiv[bestoptions]=ndiv
      ndiv
    }
    
    #number of units to be assigned to best j treatments:
    numberinbest=sapply(1:k, function(j) (min(ceiling(N/k),floor(N/j)) :floor(N/j)))
    nmatrix=sapply(1:k, function(j) sapply(numberinbest[[j]],
                                      function(nn) divider(j,nn)))
    
    unique(t(do.call(cbind, nmatrix))) #convert list to matrix with assignments in rows, remove duplicate rows
  }  
  
}




####
#for testing
# A=c(1,1,1)
# B=c(1,1,1)
# C=c(0,0,0)
# 
# Seed(A,B)
# Uhat(A,B,C,c(0,0,10))
# Uhat(A,B,C,c(10,0,0))
# Uhat(A,B,C,c(5,5,5))
# Uhat(A,B,C,c(10,10,10))
# Uhat(A,B,C,c(20,20,20))
