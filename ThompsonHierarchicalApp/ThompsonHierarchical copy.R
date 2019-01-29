betabinomialMLE = function(NN,SS) { #matrix of trials and successes, across treatments and strata
    k=dim(NN)[1]
    nx=dim(NN)[2]
    
    AB=matrix(0,k,2)
    LLH=rep(0,k)
    
    # for each treatment arm, find MLE of hyperparameters
    for (d in 1:k) {
        minusLLH=function(ab) { #negative of betabinomial log likelihood, up to constant
            -sum(mapply(lbeta, SS[d,]+ab[1], NN[d,]-SS[d,] + ab[2])) +  nx* lbeta(ab[1], ab[2])
        }
        findmin<-nlminb(objective=minusLLH, 
                        start=c(1,1),
                        lower=c(0,0),
                        upper=rep(500, 2),
                        control = list(eval.max = 10^5, 
                                       iter.max = 10^5, 
                                       abs.tol = 10^(-8)))
        AB[d,] = findmin$par
        LLH[d] = -findmin$objective
    }
    
    list(AB=AB, LLH=LLH)
}

hierarchicalPosteriorDraw = function(NN,SS, 
                                     LLH) { #vector of maximized log-likelihoods, up to constant, for each treatment arm
    k=dim(NN)[1]
    nx=dim(NN)[2]
    AB=matrix(0,k,2)
    f=rep(0,k)
    
    # could add option to the following, to simply take MLE instead of posterior draw for hyperparameter?
    # => TBD: add empirical Bayes option to function!
    
    for (d in 1:k) { #loop over treatment values
        crit=-1
        while (crit<0) { #rejection sampling from the hyper-posterior
            nn=rchisq(1,3) #chi squared prior for precision
            mn=runif(1) #uniform prior for mean
            ab=c(mn*nn, (1-mn)*nn)
            f[d]=sum(mapply(lbeta, SS[d,]+ab[1], NN[d,]-SS[d,] + ab[2])) -  k* lbeta(ab[1], ab[2])  #negative of betabinomial log likelihood, up to constant
            crit= f[d]-LLH[d]- log(runif(1))
        }
        AB[d,]=ab
    }

    list(theta=matrix(mapply(rbeta, 1,  SS+AB[,1], NN-SS + AB[,2]),k,nx), #draws from conditional beta posterior
         A=AB[,1], B=AB[,2], f=f)
}


# like posterior draw, but for posterior expectation of theta
hierarchicalPosteriorMean = function(Y,D,X,
                                     draws=1000) { #Monte Carlo draws to average
    SS=tapply(Y,list(D,X),sum) #matrix of successes
    NN=tapply(Y,list(D,X),length) #matrix of trials
    LLH=betabinomialMLE(NN,SS)$LLH
    
    k=dim(NN)[1]
    nx=dim(NN)[2]
    thetasum=matrix(0, k, nx)
    AB=matrix(0,k,2)
    f=rep(0,k)
    
    for (i in 1:draws){ 
        for (d in 1:k) { #loop over treatment values
            crit=-1
            while (crit<0) { #rejection sampling from the hyper-posterior
                nn=rchisq(1,3) #chi squared prior for precision
                mn=runif(1) #uniform prior for mean
                ab=c(mn*nn, (1-mn)*nn)
                f[d]=sum(mapply(lbeta, SS[d,]+ab[1], NN[d,]-SS[d,] + ab[2])) -  k* lbeta(ab[1], ab[2])  #negative of betabinomial log likelihood, up to constant
                crit= f[d]-LLH[d]- log(runif(1))
            }
            AB[d,]=ab
        }
        thetasum=thetasum + (SS+AB[,1]) / (NN + AB[,1] + AB[,2])
    }
    thetasum / draws
}



DtchoiceThompson=function(Y,D, #outcomes and treatments thus far
                          k, #number of treatments
                          Nt){ # number of observations for period t
  
  SS=tapply(Y,D,sum) #vector of successes
  NN=tapply(Y,D,length) #vector of trials
  A=1+SS
  B=1+NN-SS
  
  Dt=rep(0,Nt)
  previousD=-Inf # auxiliary variable to avoid repeat assignments of same D
  
  for (i  in 1:Nt) {
    thetadraw=sapply(1:k, function(j) rbeta(1, A[j], B[j]))
    Dt[i]=which.max(thetadraw)
    if (Dt[i] == previousD) {
      thetadraw[previousD] = -Inf
      Dt[i]=which.max(thetadraw)
    }
    previousD = Dt[i]
  }
  
  Dt
}



DtchoiceThompsonModified=function(Y,D, #outcomes and treatments thus far
                                  k, #number of treatments
                                  Nt, # number of observations for period t
                                  RR){ #number of replication draws
  
  Dt=rep(0,Nt)
  
  # Repeat Thompson sampling RR times
  DtRR=DtchoiceThompson(Y,D,k,Nt*RR)
  N_Dt_average=tabulate(DtRR,k) / RR #average count for each treatment value and covariate value, replicated sample
  
  
  N_Dt_floor=floor(N_Dt_average) #average number of assignments, rounded down
  Rem_D = N_Dt_average - N_Dt_floor #remainder
  n_Rem= Nt - sum(N_Dt_floor) #remaining number of units to be assigned
  
  if (n_Rem> 0) {
    Dt= c(rep(1:k, N_Dt_floor), # assigning each treatment acoording to rounded down average count
          sample(1:k,size=n_Rem,replace=T, prob=Rem_D)) #remaining units assigned randomly from remain replicate assignments
  } else {
    Dt= rep(1:k, N_Dt_floor)
  }
  
  sample(Dt) #randomly permute
}



DtchoiceThompsonHierarchical=function(Y,D,X, #outcomes, treatments, and covariates thus far
                                      k,nx, #number of treatments and number of strata
                                      Xt){ # covariates for period t
  
    SS=tapply(Y,list(D,X),sum) #matrix of successes
    NN=tapply(Y,list(D,X),length) #matrix of trials
    
    MLE=betabinomialMLE(NN,SS)
    
    Nt=length(Xt)
    Dt=rep(0,Nt)
    for (i  in 1:Nt) {
        thetadraw=hierarchicalPosteriorDraw(NN,SS,MLE$LLH)$theta[,Xt[i]] #draw from posterior for covariate value Xt[i]
        Dt[i]=which.max(thetadraw)
    }
    
  Dt
}



DtchoiceThompsonHierarchicalAlternating=function(Y,D,X, #outcomes, treatments, and covariates thus far
                                      k,nx, #number of treatments and number of strata
                                      Xt){ # covariates for period t
  
  SS=tapply(Y,list(D,X),sum) #matrix of successes
  NN=tapply(Y,list(D,X),length) #matrix of trials
  
  MLE=betabinomialMLE(NN,SS)
  
  Nt=length(Xt)
  Dt=rep(0,Nt)
  previousD=rep(-Inf, nx) # auxiliary vector to avoid repeat assignments of same D
  
  for (i  in 1:Nt) {
    thetadraw=hierarchicalPosteriorDraw(NN,SS,MLE$LLH)$theta[,Xt[i]] #draw from posterior for covariate value Xt[i]
    Dt[i]=which.max(thetadraw)
    if (Dt[i] == previousD[Xt[i]]) {
      thetadraw[Dt[i]] = -Inf
      Dt[i]=which.max(thetadraw)
    }
    previousD[Xt[i]] = Dt[i]
  }
  
  Dt
}



DtchoiceThompsonHierarchicalModified=function(Y,D,X, #outcomes, treatments, and covariates thus far
                                      k,nx, #number of treatments and number of strata
                                      Xt, # covariates for period t
                                      RR){ #number of replication draws
  
  N_Xt=table(Xt) #count of covariate cells
  Nt=length(Xt)
  Dt=rep(0,Nt)
  
  # Repeat Hierarchical Thompson sampling RR times for Xt
  XtRR=rep(Xt,RR)
  DtRR=DtchoiceThompsonHierarchicalAlternating(Y,D,X,k,nx,XtRR)
  N_XtDt_average=table(list(X=XtRR,D=DtRR)) / RR #average count for each treatment value and covariate value, replicated sample
  

  #loop over each covariate value to get assignment vectors
  for (x in 1:nx) {
    N_xDt_floor=floor(N_XtDt_average[x,]) #average number of assignments, rounded down
    Rem_xD = N_XtDt_average[x,] - N_xDt_floor #remainder
    n_Rem= N_Xt[x] - sum(N_xDt_floor) #remaining number of units to be assigned
    if (n_Rem> 0) {
      Dtx= c(rep(1:k, N_xDt_floor), # assigning each treatment acoording to rounded down average count
            sample(1:k,size=n_Rem,replace=T, prob=Rem_xD)) #remaining units assigned randomly from remain replicate assignments
    } else {
      Dtx= rep(1:k, N_xDt_floor)
    }  
    Dt[Xt==x] = sample(Dtx) #randomly permute
  }
  
  
  Dt
}




SimulateY=function(theta, D, X){
    N=length(X)
    thetaDX=mapply(function(d,x) theta[d,x], D,X)
    Y=runif(N)<thetaDX    
}

SimulateX=function(PX,N){
  nx=length(PX)
  sample(x = 1:nx, N, replace = TRUE, prob = PX) 
}



#divide units for each value of X equally across treatments
StratifiedAssignment=function(X,k,nx){
  N=length(X)
  D=rep(0,N)
  nextD=sample(1:k,nx, replace = TRUE) #random starting values in each stratum
  for (i in 1:N) {
    D[i]=nextD[X[i]]
    nextD[X[i]]= (nextD[X[i]] %% k ) + 1 #rotating through treatment values
  }
  D
}



DtchoiceCovariates=function(Y,D,X, #outcomes, treatments, and covariates thus far
                            k,nx, #number of treatments and number of strata
                            Xt, # covariates for period t
                            method="stratified"){
  if (method=="stratified") {
    Dt=StratifiedAssignment(Xt,k,nx)
  } else if (method=="random") {
    Dt=StratifiedAssignment(X=rep(0,length(Xt)),1,nx)
  } else if (method=="thompson") {
    Dt=DtchoiceThompsonHierarchical(Y,D,k,Nt)
  } else if (method=="modifiedthompson") {
    Dt=DtchoiceThompsonHierarchicalModified(Y,D,k,Nt, RR=5)
  } 
  
  Dt
}  


