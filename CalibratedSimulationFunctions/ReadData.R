DataToTheta=function(filename, dataname, k, strataVars, outcomename,treatmentname, covariatesnames, printFigures=FALSE){
    Data=read_csv(paste("../Datasets/Cleaned/", filename, ".csv", sep=""))
    N=nrow(Data)
    
    #check for missings?
    
    treatDummies=paste("treatment",1:k, sep="") #names of treatment variables
    
    Data=Data %>%
        mutate(treatment=factor(as.matrix(Data[treatDummies])%*%(1:k))) %>%
        mutate(Strata=(interaction(select(.,strataVars))))
    
    #recoding the levels in a reproducible way    
    oldlevels=levels(Data$Strata)  
    key=data_frame(Strata=oldlevels, strata=factor(1:length(oldlevels)))
    Data=Data %>%
        left_join(., key, by = "Strata") %>%
        select(-Strata)
    
    #average outcomes by treatment and stratum
    sumstats=Data %>% 
        group_by(treatment, strata) %>%
        summarize(meanout=mean(outcome), sumout=sum(outcome), obs=n()) 
    
    #average outcomes by treatment alone
    theta=Data %>% 
        group_by(treatment) %>%
        summarize(meanout=mean(outcome), obs=n())
    
    stratasizes = sumstats %>%
        group_by(strata) %>%
        summarize(n=sum(obs))
    
    
    # Figures and tables
    if (printFigures)
        PrintDataFigures(stratasizes,sumstats,theta, filename,dataname,outcomename,treatmentname, k)
    
    
    
    # converting to wide matrices of successes and trials, dropping NA strata
    nstrata=nrow(key)
    
    SS= sumstats  %>% #number of successes, treatments by row, strata by column
        select(treatment, strata, sumout) %>%
        spread(strata, sumout) %>%
        ungroup %>% 
        select(paste(1:nstrata)) %>%
        data.matrix()
    
    NN= sumstats  %>% #number of trials, treatments by row, strata by column
        select(treatment, strata, obs) %>%
        spread(strata, obs) %>%
        ungroup %>% 
        select(paste(1:nstrata)) %>%
        data.matrix()
    
    NX= stratasizes %>% #number of units per stratum
        slice(1:nstrata) %>%
        select(n) %>%
        data.matrix()
    
    PX=NX/sum(NX)
    
    
    list(filename=filename, dataname=dataname, N=N,
         theta=theta$meanout, sumstats=sumstats, stratasizes =stratasizes, key=key, #output for old simulations (no covariates)
         SS=SS, NN=NN, PX=PX, k=k) #output for thompson simulations (using covariates)
}



#produce figures and tables of calibrated parameter values
PrintDataFigures=function(stratasizes,sumstats,theta, filename,dataname,outcomename,treatmentname, k){
    

    FooterText=paste("Outcome: ", outcomename, 
                   "\nTreatments: ", treatmentname, 
                   ".", sep="")
  
  
    #careful: we are dropping "missing" strata from figures!
    ggplot(drop_na(stratasizes),aes(x=factor(strata, levels = rev(levels(strata))), y=n)) +
        geom_col(fill=fillcolor, width=.2) + 
        theme_light() +
        theme(#panel.background = element_rect(fill = backcolor, colour = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
        coord_flip() +
        labs(title=dataname,
             x="Stratum", y="Number of observations")
    
    ggsave(paste("../Figures/Applications/", filename, "strata.pdf", sep=""), width = 4, height =1+ .3*length(levels(sumstats$strata)))
    
    xmax=1
    if (max(sumstats$meanout) < .5) xmax= max(sumstats$meanout)+.02
    

    
    ggplot(drop_na(sumstats), aes(x=meanout, y=factor(treatment,levels=c(k:1)))) +
        geom_point(color=fillcolor) +
        #geom_segment(aes(x=0, xend=meanout, yend=factor(treatment,levels=c(k:1))), color=fillcolor, size=.5) +
        #scale_y_discrete(labels=paste("treatment", 1:k)) +
        facet_grid(strata~.) +
        scale_x_continuous(limits=c(0,xmax)) +
        theme_light() +
        theme(#panel.background = element_rect(fill = backcolor, colour = NA),
            #panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.x=element_blank(),
            plot.caption=element_text(size=7)) +
        labs(title=dataname,
             x="Mean outcome", y="Treatment",
             caption=FooterText)
    
    ggsave(paste("../Figures/Applications/", filename, ".pdf", sep=""), width = 4, height = 0.15*length(levels(sumstats$strata))*k+1.5)

    

    xmax=1
    if (max(theta$meanout) < .3) xmax= max(theta$meanout)+.02
    
    ggplot(theta, aes(x=meanout, y=factor(treatment,levels=c(k:1)))) +
      geom_point(color=fillcolor) +
      #geom_segment(aes(x=0, xend=meanout, yend=factor(treatment,levels=c(k:1))), color=fillcolor, size=.5) +
      #scale_y_discrete(labels=paste("treatment", 1:k)) +
      scale_x_continuous(limits=c(0,xmax)) +
      theme_light() +
      theme(#panel.background = element_rect(fill = backcolor, colour = NA),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        plot.caption=element_text(size=7)) +
      labs(title=dataname,
           x="Mean outcome", y="Treatment",
           caption=FooterText)
    
    ggsave(paste("../Figures/Applications/", filename, "_NoStrata.pdf", sep=""), width = 4, height = 0.15*k+1.5)
    
        
    write_csv(sumstats, paste("../Figures/Applications/", filename, "Sumstats.csv", sep=""))
    
}

# Alternative, more compact printout of parameter values
print_one_datafigure=function(DataList) {
    thetas=map(DataList, "theta")
    names=map(DataList, "dataname")
    NoTreatments=map(thetas,length)
    
    
    plot_data=tibble(thetas = unlist(thetas),
                     names = unlist(rep(names, NoTreatments)))
                     
        
    ggplot(plot_data, aes(x=fct_rev(as.factor(names)), y=thetas)) +
        geom_point(color=fillcolor,size=2,alpha=.7) +
        scale_y_continuous(limits=c(0,1)) +
        coord_flip() +
        theme_light() +
        theme(#panel.background = element_rect(fill = backcolor, colour = NA),
            #panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.x=element_blank(),
            plot.caption=element_text(size=7)) +
        labs(#title="Calibrated parameter values",
             x="", y="Average outcome for each treatment")
    
    ggsave("../Figures/Applications/CalibratedTheta_NoStrata.pdf", width = 8, height=2.5)    
    
    for (application in 1:length(DataList)) {
        
        ggplot(DataList[[application]]$sumstats, 
               aes(x=factor(strata, levels = rev(levels(strata))), y=meanout,group=treatment, colour=treatment)) +
            geom_line(size=.2,alpha=.5) +
            geom_point(size=2,alpha=.7) +
            scale_y_continuous(limits=c(0,1)) +
            scale_colour_viridis_d() +
            guides(colour=FALSE) +
            coord_flip() +
            theme_light() +
            theme(#panel.background = element_rect(fill = backcolor, colour = NA),
                #panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.ticks.x=element_blank(),
                plot.caption=element_text(size=7)) +
            labs(title=DataList[[application]]$dataname,
                 subtitle="Calibrated parameter values",
                 x="Stratum", y="Mean outcome for each treatment")
        
        
        ggsave(paste("../Figures/Applications/Compact_", DataList[[application]]$filename, ".pdf", sep=""), 
               width = 8, height = 0.3*length(DataList[[application]]$PX) +1)
    }
}


ReadAllData=function(printFigures=F){
    #read table of all applications from csv file
    ApplicationTable=read_csv("../Datasets/ApplicationTable.csv")
  
    #read in data for each row of Application table, store output in each element of DataList
    DataList=pmap(ApplicationTable, function(filename, dataname, k, stratavars, outcomename,treatmentname, covariatesnames) #make sure ApplicationTable has these column names
                  DataToTheta(filename, dataname, k, stratavars, outcomename, treatmentname, covariatesnames, printFigures=printFigures))
        
    if (printFigures) print_one_datafigure(DataList)
    
    DataList
}
