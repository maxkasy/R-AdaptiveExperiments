#This script implements calibrated simulations, sampling outcomes from the empirical distribution of an actual experiment.
library(future.apply)
library(tidyverse)
library(forcats)
library(xtable)
library(magrittr)

set.seed(12231983)
# prior precision
alpha0=1



SimulateTWaveDesign=function(wavesizes,C,theta, method="modifiedthompson"){
  k=length(theta) #number of treatment arms
  theta=sample(theta) #randomly permute so as not to privilege any options in policy choice

  if (method=="conventional") {
      D=EqualAssignment(sum(wavesizes),k)
      Y=simulatedSample(D,theta) #bernoulli draws with probability theta(D) - drawing from sample distribution
  } else {
      waves=length(wavesizes) #number of waves
      MM=cumsum(wavesizes)
      
      D=rep(0,MM[waves])
      Y=rep(0,MM[waves])
      
      Dt=EqualAssignment(wavesizes[1],k)
      D[1:wavesizes[1]]=Dt
      Y[1:wavesizes[1]]=simulatedSample(Dt,theta) #bernoulli draws with probability theta(D) - drawing from sample distribution
  
      for (t in 2:waves) {
        posterior=betaposterior(D[1:MM[t-1]],Y[1:MM[t-1]])
        Dt=Dtchoice(posterior$A,posterior$B,C, wavesizes[t], method)
        D[(MM[t-1]+1):MM[t]]=Dt
        Y[(MM[t-1]+1):MM[t]]=simulatedSample(Dt,theta) #bernoulli draws with probability theta(D) - drawing from sample distribution
      }
  }
  Regret(D,Y,C,theta)
}


ExpectedRegret=function(wavesizes,C,theta,methods,R){
  k=length(theta) #number of treatment arms
  
  #parallelize simulations
  plan(multiprocess)
  
  Methods=c("conventional",
                   "optimalhandpicked",
                   "optimalrandom",
                   "optimalhat",
                   "optimalhatrandom",
                   "optimal",
                   "besthalf",
                   "thompson",
                   "expectedthompson",    
                   "modifiedthompson",
                   "toptwo",
                   "modifiedthompsonanalytic",
                   "weightedmodifiedthompson")   
  MethodNames=c("non-adaptive",
                 "handpicked",
                 "random",
                 "estimated U",
                 "estimated U, random set",
                 "optimized",
                 "best half",
                 "Thompson",
                 "expected Thompson",
                 "alternating Thompson",
                 "top-two Thompson",
                 "modified Thompson",
                 "weighted modified Thompson")
  shareTreatments=list() #empty list to store vectors of shares assigned to each treatment
      
  for (i in methods) { #pick here which methods to simulate
      sink("status_ExpectedRegret.txt") #status file to track computations
      cat("waves ", wavesizes, "\n", 
          "theta", theta, "\n",
          "method", Methods[i])
      sink()  
    
    regretTWave=future_sapply(1:R, function(j) SimulateTWaveDesign(wavesizes,C,theta, Methods[i]))
    regretTable=rbind(get0("regretTable"),
                      tibble(Statistic=paste("$\\quad$ ", MethodNames[i], sep=""),
                                Value=mean(regretTWave[1,])))
    shareoptTable=rbind(get0("shareoptTable"),
                        tibble(Statistic=paste("$\\quad$ ", MethodNames[i], sep=""),
                                Value=mean(regretTWave[1,]==0)))
    
    shareTreatments[[Methods[i]]]=table(factor(regretTWave[1,], levels=max(theta)-theta)) #store shares assigned to each treatment for method i
  }
  
  lastrows=tibble(Statistic=c("Units per wave"),#, "Number of treatments"),
                 Value=c(wavesizes[1]))#, k))
  
  tablecolumns=rbind(tibble(Statistic="Regret", Value=NA),
                     regretTable, 
                     tibble(Statistic="Share optimal", Value=NA),
                     shareoptTable, 
                     lastrows)

  list(tablecolumns = tablecolumns, shareTreatments=shareTreatments)
}



DesignTable=function(DataList,methods,MC_replicates=100,columnames=NULL,filename=NULL, DoHistograms=F) {

  # Run Expected Regret simulations for each value of theta    
  ResultsTemp=map(DataList, function(data) 
                    ExpectedRegret(data$wavesizes,C=rep(0, length(data$theta)),data$theta,methods,MC_replicates))      
    
  RegretTableTemp=map(ResultsTemp, "tablecolumns") #extract table columns
  shareTreatmentsList=map(ResultsTemp, "shareTreatments") #extract shares assigned to treatments
  
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
 
  #write to tex file, plot histograms
  if (!is.null(filename)) {
      PrintRegretTable(RegretTable,filename,DataList[[1]]$dataname,MC_replicates, length(methods))
      if (DoHistograms) PrintRegretHistogram(shareTreatmentsList,filename, DataList[[1]]$dataname, #paste(DataList[[1]]$filename, " et al.", sep=""), 
                                             MC_replicates, map_dbl(DataList, "waves"))
  }      
}



PrintRegretTable = function(RegretTable,filename,caption, MC_replicates, nmethods){
    #write_csv(RegretTable, paste("../Figures/", filename, "RegretTable.csv", sep=""))
    labtext=paste(filename,"_", 
                  #paste(wavesizes, collapse="_"),"_",
                  "RegretTable", sep="")
    rows=2*(nmethods+1)
    remainderrows=dim(RegretTable)[1]-rows
    cols=dim(RegretTable)[2]
    digs=matrix(c(rep(3, rows*(cols+1)), rep(0,remainderrows*(cols+1))), nrow=rows+remainderrows, ncol=cols+1, byrow=T) #controling digits
    
    print(xtable(RegretTable, type = "latex",
                 caption=caption,#paste(MC_replicates, " replications, ", waves, " waves.", sep=""),
                 label=paste("tab:", labtext, sep=""),
                 digits=digs),
          hline.after = c(-1,0,nmethods+1,rows,rows+remainderrows), #horizontal lines
          file = paste("../Figures/Simulations/", labtext,".tex", sep=""),
          caption.placement = "top",
          latex.environments = "center", #centering the table and caption
          include.rownames=FALSE, #to omit row numbering
          sanitize.text.function=function(x){x}) #to maintain tex strings as intended
}




PrintRegretHistogram=function(shareTreatmentsList,filename, dataname,MC_replicates, waves) {
  histTibble=tibble()
  for (i in 1:length(columnames)) {
    shareTreatments=shareTreatmentsList[[i]]
    
    tibble_i=tibble(regret=as.double(names(shareTreatments[[1]])),
                    non_adaptive= as.double(shareTreatments[[length(shareTreatments)]] / MC_replicates),
                    modifiedThompson= as.double(shareTreatments[[1]] / MC_replicates),
                    waves=columnames[i])
    #redundant row to make CDFs come out nicely
    tibble_i %<>% add_row(regret=-0.001, non_adaptive=0, modifiedThompson=0, waves=columnames[i])    
    
    histTibble %<>% bind_rows(
      gather(data=tibble_i,
             key=algorithm,
             value=probability,
             c(modifiedThompson,non_adaptive) )
    )
  }
  histTibble$waves %<>% factor(levels=c("2 waves", "4 waves", "10 waves"))
  histTibble$algorithm %<>% factor(levels=c("non_adaptive", "modifiedThompson"))
  maxregret=max(histTibble$regret)

  #calculate CDF of regret
  histTibble %<>% arrange(waves, algorithm,regret, probability) %>%
    group_by(waves, algorithm) %>%
    mutate(cdf=cumsum(probability))
  
  #print histograms
  ggplot(histTibble[histTibble$regret>=0,], aes(x=regret, y=probability, fill=algorithm)) +
    geom_col(width=maxregret/15, 
             position=position_dodge2(padding=.2, preserve = c("single"))) +
    scale_fill_manual(labels = c("non-adaptive", "modified Thompson"),
                      values = c("skyblue4", "skyblue")) +
    scale_y_continuous(limits=c(0,1)) +
    #scale_x_reverse() + #(breaks=unique(histTibble$regret)) +
    facet_grid(cols=vars(waves) ) +
    #scale_fill_viridis(discrete=TRUE) +
    coord_flip() +
    theme_light() +
    theme(#panel.background = element_rect(fill = backcolor, colour = NA),
      #panel.grid.major = element_blank(),
      #panel.grid.minor = element_blank(),
      axis.ticks.x=element_blank(),
      axis.ticks.y=element_blank(),
      #plot.caption=element_text(size=7),
      legend.position = "top",
      strip.text = element_text(size = 10)) +
    labs(title=dataname, fill="",
         x="Regret", y="Share of simulations")
  
  hist_pathname=paste("../Figures/Simulations/Histograms/",
                      filename, "_RegretHistogram.pdf", sep="")
  ggsave(hist_pathname, width = 8, height = 3.5)
  
  
  # print quantiles
  ggplot(histTibble, aes(x=regret, y=cdf, color=algorithm)) +
    geom_step(size=1) +
    scale_color_manual(labels = c("non-adaptive", "modified Thompson"),
                       values = c("skyblue4", "skyblue")) +
    scale_y_continuous(limits=c(0,1)) +
    facet_grid(cols=vars(waves) ) +
    #scale_fill_viridis(discrete=TRUE) +
    coord_flip() +
    theme_light() +
    theme(axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          #plot.caption=element_text(size=7),
          legend.position = "top",
          strip.text = element_text(size = 10)) +
    labs(title=" ", 
         color="",
         x="Quantile of regret", 
         y="Share of simulations")    
  
  quantile_pathname=paste("../Figures/Simulations/Histograms/",
                          filename, "_RegretQuantile.pdf", sep="")
  ggsave(quantile_pathname, width = 8, height = 3.5)
}



# 
# PrintRegretHistogram=function(shareTreatmentsList,filename, dataname,MC_replicates, waves) {
#     for (i in 1:length(columnames)) {
#         pathname=paste("../Figures/Simulations/Histograms/",
#                        filename,"_", gsub(" ", "", columnames[i]),
#                        "_RegretHistogram.pdf", sep="")
#         
#         shareTreatments=shareTreatmentsList[[i]]
#         
#         histTibble=tibble(regret=as.double(names(shareTreatments[[1]])),
#                            non_adaptive= as.double(shareTreatments[[length(shareTreatments)]] / MC_replicates),
#                            modifiedThompson= as.double(shareTreatments[[1]] / MC_replicates)) #careful about coordinates being right!
# 
#         ggplot(histTibble, aes(x=regret, y=modifiedThompson)) +
#             geom_point(color=fillcolor, size=1) +
#             geom_segment(aes(x=regret,xend=regret, y=non_adaptive,  yend=modifiedThompson), color=fillcolor, size=.5) +
#             scale_y_continuous(limits=c(0,1)) +
#             #coord_cartesian(ylim=c(0,1)) +
#             theme_light() +
#             theme(#panel.background = element_rect(fill = backcolor, colour = NA),
#                 #panel.grid.major = element_blank(),
#                 panel.grid.minor = element_blank(),
#                 axis.ticks.x=element_blank(),
#                 axis.ticks.y=element_blank(),
#                 plot.caption=element_text(size=7)) +
#             labs(title=paste(dataname, ", ", waves[i], " waves", sep=""),
#                  x="Regret", y="Share of simulations")
#         
# 
#         
#         ggsave(pathname, width = 4, height = 4)
#     }
# }
# 



