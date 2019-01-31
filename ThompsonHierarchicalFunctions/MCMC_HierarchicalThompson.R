################################
# Thompson with covariates
# helper functions for MCMC

# prior for hyperparameters governing distribution of theta across strata within each treatment arm 
log.prior = function(alpha,beta) {
  -2.5*log(max(alpha + beta,0))
}

# sampling theta vector from posterior, given hyperparameters, for a given treatment arm
draw.thetas = function(alpha,beta, NNd, SSd, nx) {
  rbeta(nx,alpha+SSd,beta+NNd-SSd)
}

# sampling from posterior for hyperparameters, given theta vector, for a given treatment arm 
draw.alpha = function(alpha,beta,theta,prop.sd,nx) {
  alpha.star = max(rnorm(1,alpha,prop.sd),.5)
  num = nx*(lgamma(alpha.star+beta) - lgamma(alpha.star)) +
    alpha.star*sum(log(theta)) + log.prior(alpha.star,beta)
  den = nx*(lgamma(alpha+beta)      - lgamma(alpha)) +
    alpha     *sum(log(theta)) + log.prior(alpha,beta)
  acc = ifelse((log(runif(1))<=num - den)&&(alpha.star>0),1,0)

  ifelse(acc,alpha.star,alpha)
}

draw.beta = function(alpha,beta,theta,prop.sd,nx) {
  beta.star = max(rnorm(1,beta,prop.sd),.5)
  num = nx*(lgamma(alpha+beta.star) - lgamma(beta.star)) +
    beta.star*sum(log(1-theta)) + log.prior(alpha,beta.star)
  den = nx*(lgamma(alpha+beta)      - lgamma(beta)) +
    beta     *sum(log(1-theta)) + log.prior(alpha,beta)
  acc = ifelse((log(runif(1))<=num - den)&&(beta.star>0),1,0)
  
  ifelse(acc,beta.star,beta)
}

sample.theta.d = function(NNd, SSd, nx, 
                          RR=2000) { #sampling period
  B = 1000 #burn in period
  MM = B + RR
  # Metropolis tuning parameters
  alpha.prop.sd =  1
  beta.prop.sd =   1
  
  alpha = rep(0,MM)
  beta = alpha
  theta = matrix(0,MM,nx)
  
  # Initial values for the chain
  alpha[1] = 2
  beta[1] = 2
  theta[1,] = draw.thetas(alpha[1],beta[1], NNd, SSd,nx)
 
  # MCMC simulation
  for (m in 2:MM) {
    alpha[m] = draw.alpha(alpha[m-1],beta[m-1],theta[m-1,],alpha.prop.sd,nx)
    beta[m] = draw.beta(alpha[m],beta[m-1],theta[m-1,],beta.prop.sd,nx)
    theta[m,] = draw.thetas(alpha[m],beta[m], NNd, SSd,nx)
  }

  theta[(B+1):MM,]
}

DtchoiceMCMCProbabilities=function(Y,D,X, #outcomes, treatments, and covariates thus far
                                   k,nx, #number of treatments and number of strata
                                   C=rep(0,k), #vector of treatment cost
                                   RR=2000){ #number of replication draws
  
  SS=tapply(Y,list(D,X),sum, default=0) #matrix of successes
  NN=tapply(Y,list(D,X),length, default=0) #matrix of trials

  P_Dt=matrix(0,nx,k)
  thetadraws=list()

  for (d in 1:k) {
    thetadraws[[d]]=sample.theta.d(NN[d,], SS[d,], nx, RR)
  }
  
  Dt_x=factor(rep(0,RR), levels=1:k)
  thetaxdraw=rep(0,k)
  for (x in 1:nx) {
    for (r in 1:RR) {
      for (d in 1:k) thetaxdraw[d]=thetadraws[[d]][r,x]
      Dt_x[r]=which.max(thetaxdraw-C)
    }
    P_Dt[x,]=tabulate(Dt_x, nbins=k)/RR
  }
  
  P_Dt=as_tibble(P_Dt)
  colnames(P_Dt)=paste(1:k)
  P_Dt
}



hierarchicalPosteriorMean = function(Y,D,X,k,nx) {
    SS=tapply(Y,list(D,X),sum, default=0) #matrix of successes
    NN=tapply(Y,list(D,X),length, default=0) #matrix of trials

    thetadraws=list()
    
    for (d in 1:k) {
        thetadraws[[d]]=sample.theta.d(NN[d,], SS[d,], nx)
    }
    
    thetameans=matrix(0,k,nx)
    
    for (d in 1:k) {
        thetameans[d,]=colMeans(thetadraws[[d]])
    }
    thetameans
}


##########################################################
# Treatment assignment functions

#divide units for each value of X equally across treatments
StratifiedAssignment=function(X,k,nx){
    N=length(X)
    D=rep(0,N)
    nextD=sample(1:k,nx, replace = TRUE) #random starting values in each stratum
    for (i in 1:N) {
        D[i]=nextD[X[i]]
        nextD[X[i]]= (nextD[X[i]] %% k ) + 1 #rotating through treatment values
    }
    factor(D,levels=1:k)
}

#D vector based on shares
ProportionalAssignment = function(Shares, nt) {
    k=length(Shares)
    n_floor=floor(nt*Shares) #number of units assigned to each treatment, rounded down
    remainder = nt*Shares - n_floor
    units_remaining= nt-sum(n_floor) #remaining number of units to be assigned
    
    if (units_remaining > 0) {
        Dt= c(rep(1:k, n_floor), # assigning each treatment acoording to rounded down average count
              sample(1:k,size=units_remaining,replace=T, prob=remainder)) #remaining units assigned randomly from remain replicate assignments
    } else {
        Dt= rep(1:k, n_floor)
    }
    sample(Dt)
}


DtchoiceThompsonHierarchical=function(Y,D,X, #outcomes, treatments, and covariates thus far
                                      k,nx, #number of treatments and number of strata
                                      Xt){ # covariates for period t
    
    P_Dt=DtchoiceMCMCProbabilities(Y,D,X,k,nx)
    
    Nt=length(Xt)
    Dt=rep(0,Nt)
    for (i  in 1:Nt) {
        Dt[i]=sample(1:k,size=1, prob=P_Dt[Xt[i],])
    }
    
    factor(Dt, levels=1:k)
}


DtchoiceThompsonHierarchicalExpected=function(Y,D,X, #outcomes, treatments, and covariates thus far
                                      k,nx, #number of treatments and number of strata
                                      Xt){ # covariates for period t
    
    P_Dt=DtchoiceMCMCProbabilities(Y,D,X,k,nx) #thompson probabilities
    Xt_counts=tabulate(Xt, nbins=nx) #count of units in each stratum in period t
    
    Nt=length(Xt)
    Dt=rep(0,Nt)
    
    Dt_list=map(1:nx, #for each stratum...
                function(x) ProportionalAssignment(as_vector(P_Dt[x,]), Xt_counts[x])) #do proportional assignment based on row of P_Dt

    #spread out the proportional assignments based on covariates
    stratacounts=rep(1,nx)
    for (i  in 1:Nt) {
        Dt[i]=Dt_list[[Xt[i]]][stratacounts[Xt[i]]]
        stratacounts[Xt[i]]=stratacounts[Xt[i]]+1
    }
    
    factor(Dt, levels=1:k)
}



DtchoiceThompsonHierarchicalModified=function(Y,D,X, #outcomes, treatments, and covariates thus far
                                              k,nx, #number of treatments and number of strata
                                              Xt){ # covariates for period t
    
    P_Dt=DtchoiceMCMCProbabilities(Y,D,X,k,nx) #thompson probabilities
    Xt_counts=tabulate(Xt, nbins=nx) #count of units in each stratum in period t
    
    Nt=length(Xt)
    Dt=rep(0,Nt)
    
    Dt_list=map(1:nx, #for each stratum...
                function(x) {
                    P_Dt_x_modified=P_Dt[x,]*(1-P_Dt[x,]) #calculate modified Thompson probabilities
                    if (sum(P_Dt_x_modified) >0) P_Dt_x_modified=P_Dt_x_modified/sum(P_Dt_x_modified)
                        else P_Dt_x_modified=P_Dt[x,]
                    ProportionalAssignment(P_Dt_x_modified, Xt_counts[x]) #do proportional assignment based on modified probabilities
                })

    #spread out the proportional assignments based on covariates
    stratacounts=rep(1,nx)
    for (i  in 1:Nt) {
        Dt[i]=Dt_list[[Xt[i]]][stratacounts[Xt[i]]]
        stratacounts[Xt[i]]=stratacounts[Xt[i]]+1
    }
    
    factor(Dt, levels=1:k)
}




DtchoiceCovariates=function(Y,D,X, #outcomes, treatments, and covariates thus far
                            k,nx, #number of treatments and number of strata
                            Xt, # covariates for period t
                            method="stratified"){
    if (method=="stratified") {
        Dt=StratifiedAssignment(Xt,k,nx)
    } else if (method=="random") {
        Dt=StratifiedAssignment(X=rep(1,length(Xt)),k,1)
    } else if (method=="thompson") {
        Dt=DtchoiceThompsonHierarchical(Y,D,X, k,nx,Xt)
    } else if (method=="thompsonexpected") {
        Dt=DtchoiceThompsonHierarchicalExpected(Y,D,X, k,nx,Xt)
    } else if (method=="modifiedthompson") {
        Dt=DtchoiceThompsonHierarchicalModified(Y,D,X, k,nx,Xt)
    }
    
    Dt
}  

