SimulateY=function(theta, D, X){
    N=length(X)
    thetaDX=mapply(function(d,x) theta[d,x], D,X)
    Y=runif(N)<thetaDX    
}

SimulateX=function(PX,N){
    nx=length(PX)
    factor(sample(nx, N, replace = TRUE, prob = PX),levels=1:nx)
}

SimulateTWaveDesignCovariates=function(wavesizes,C,theta, PX, method="stratified"){
    k=length(C) #number of treatment arms
    nx = length(PX)
    theta=theta[sample(k),] #randomly permute rows so as not to privilege any options in policy choice
    
    waves=length(wavesizes) #number of waves
    MM=cumsum(wavesizes)
    
    X=SimulateX(PX,MM[waves])  #random sampling from distribution PX
    #X=1:MM[waves] %% nx + 1 #making life easier by "stratified sampling" from X
   
    D=rep(0,MM[waves])
    Y=rep(0,MM[waves])
    
    if (method=="random") {
        Dt=StratifiedAssignment(rep(1,wavesizes[1]), 1, nx)
    } else{
        Dt=StratifiedAssignment(X[1:wavesizes[1]], k, nx)
    }    
    Xt=X[1:wavesizes[1]]
    D[1:wavesizes[1]]=Dt
    Y[1:wavesizes[1]]=SimulateY(theta, Dt, Xt)
    
    for (t in seq(2, length=max(0,waves-1))) {
        previous=1:MM[t-1]
        current=(MM[t-1]+1):MM[t]
        Dt=DtchoiceCovariates(Y[previous], D[previous], X[previous],
                              k, nx, X[current], method)
        D[current]=Dt
        Y[current]=SimulateY(theta, Dt, X[current]) 
    }
    
    thetahat=hierarchicalPosteriorMean(Y,D,X,k,nx)
    
    #careful here - giving priority to lower indices...
    Dstar=apply(thetahat,2, which.max)
    regretX=apply(theta,2, max) - theta[cbind(Dstar, 1:nx)]
    regret=sum(regretX*PX)
    
    #list(X=X, D=D, Y=Y, thetahat=thetahat, Dstar=Dstar, regretX=regretX, regret=regret)
    regret
}




ExpectedRegretCovariates=function(wavesizes,C,theta,
                                  PX,methods,R){
  k=nrow(theta) #number of treatment arms
  nx=ncol(theta)
  
  Methods=c("random",
            "stratified",
            "thompson",
            "thompsonexpected",
            "modifiedthompson")   
  MethodNames=c("non-adaptive",
                "non-adaptive stratified",
                "Thompson",
                "expected Thompson",
                "modified Thompson")
  shareTreatments=list() #empty list to store vectors of shares assigned to each treatment
#browser()
  #parallelize simulations
  no_cores = detectCores()
  clust = makeCluster(no_cores, type="FORK")  #forking requires mac or Linux OS!
    
  for (i in methods) { #pick here which methods to simulate
    sink("status_ExpectedRegret_Covariates.txt") #status file to track computations
    cat("waves ", wavesizes, "\n",
       "theta", theta, "\n",
       "method", Methods[i])
    sink()

    regretTWave=parSapply(clust, 1:R, function(j) SimulateTWaveDesignCovariates(wavesizes,C,theta,PX, Methods[i]))
    regretTable=rbind(get0("regretTable"),
                      tibble(Statistic=paste("$\\quad$ ", MethodNames[i], sep=""),
                             Value=mean(regretTWave)))
    shareoptTable=rbind(get0("shareoptTable"),
                        tibble(Statistic=paste("$\\quad$ ", MethodNames[i], sep=""),
                               Value=mean(regretTWave==0)))
    #shareTreatments[[Methods[i]]]=table(factor(regretTWave)) #store shares assigned to each treatment for method i
  }
  
  
  stopCluster(clust)
  
  lastrows=tibble(Statistic=c("Units per wave"),Value=c(wavesizes[1]))
  # lastrows=tibble(Statistic=c("Units per wave", "Number of treatments", "Number of strata"),
  #                 Value=c(wavesizes[1], k, nx))

  tablecolumns=rbind(tibble(Statistic="Regret", Value=NA),
                     regretTable, 
                     tibble(Statistic="Share optimal", Value=NA),
                     shareoptTable, 
                     lastrows)
  
  list(tablecolumns = tablecolumns)#, shareTreatments=shareTreatments)
}




DesignTableCovariates=function(DataList,methods,MC_replicates=100,columnames=NULL,filename=NULL) {
  # Run Expected Regret simulations for each value of theta    
  ResultsTemp=map(DataList, function(data) 
    ExpectedRegretCovariates(data$wavesizes,C=rep(0, length(data$theta)),theta=data$SS/data$NN,data$PX,methods,MC_replicates))      
  
  RegretTableTemp=map(ResultsTemp, "tablecolumns") #extract table columns
  #shareTreatmentsList=map(ResultsTemp, "shareTreatments") #extract shares assigned to treatments
  
  # Combine tables
  RegretTable=bind_cols(RegretTableTemp[[1]]["Statistic"], #row labels
                        map(RegretTableTemp, 2)) #keep only simulation values
  #change column names
  if (!is.null(columnames)) {
    colnames(RegretTable)=c("Statistic", columnames)
  } else {
    colnames(RegretTable)=c("Statistic", map(DataList, function(data) 
      paste("$\\theta = (", toString(data$theta,sep=", "), ")$")))
  }
  
  #write to tex file
  waves=length(DataList[[1]]$wavesizes)
  if (!is.null(filename)) {
    PrintRegretTable(RegretTable,filename, DataList[[1]]$dataname,MC_replicates, length(methods))
  }      
}





# RunAllSimulationsThompson=function(DataList,
#                                    waves = 4, #number of waves
#                                    nt = 36, #units per wave
#                                    RR = 1000){ #number of replications for simulation
# 
#     simRegret=matrix(0,2,2) #columns for applications, rows for methods
#     
#     for (application in 1:2){
#         theta=DataList[[application]]$SS / DataList[[application]]$NN
#         PX=DataList[[application]]$PX
#         #note: need enough units to observe each treatment/covariate combo in first round, for MLE
# 
#         C=rep(0,DataList[[application]]$k)
#         
#         for (r in 1:RR) {
#             #note 2: permute theta rows for fair comparisons
#             
#             # Thompson design
#             simDesign=SimulateTWaveDesignThompson(wavesizes = rep(nt,waves),C,theta, PX)
#             # 1 wave comparison with the same number of units. stratified assignment.
#             simDesign1wave=SimulateTWaveDesignThompson(wavesizes = nt * waves,C,theta, PX)
#             
#             simRegret[1,application] = simRegret[1,application] +  simDesign$regretX %*% PX
#             simRegret[2,application] = simRegret[2,application] +  simDesign1wave$regretX %*% PX
#             
#         }    
#     }
#     simRegret=simRegret / RR
#     print(simRegret)
#     ## TBD: formatting, design choices, further comparisons?
# }



#modify this to get repeated draws, and adjust X distribution / weighting of regret
# include stratified equalsplit as comparison, and fully random assignment.

