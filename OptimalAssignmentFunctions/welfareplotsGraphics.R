#library(tidyverse)

backcolor="azure2" #background color for plots

PlotSimplex=function(A,B,C,N){
  # can choose Uhat or U for Ufunction. If Uhat, need to Seed.
  # Seed(A,B)
  USimplex=UoverSimplex(A,B,C,N, U)
  
  Uround=signif(USimplex$U+.0000001,digits=4)
  USimplex$maximizer=(Uround==max(Uround))
  USimplex$fontf=sapply(USimplex$maximizer, function(l) (if (l) "bold" else "plain"))

  plotTitle=bquote(alpha ~" = (" ~ 
                       .(paste(A, collapse=", ")) ~ 
                       "), "  ~ beta ~ " = ("  ~
                       .(paste(B, collapse=", ")) ~ 
                       ")")      
  ggplot(USimplex, aes(x=as.integer(n1), y=as.integer(n2), z=U)) +
    geom_tile(aes(fill = U))  + 
    scale_fill_gradient( low = "skyblue4", high = "white") +
    geom_text(aes(label=format(U+.0000001, digits=3), #note addition of small number is to generate consistent rounding
                  fontface=fontf),
              size=min(15/N,4),
              show.legend=FALSE) +
    coord_fixed() +
    xlab(expression(n[1])) + ylab(expression(n[2])) +
    scale_x_continuous(breaks = 0:N) + scale_y_continuous(breaks = 0:N) +
    theme(panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = backcolor, colour = NA),
          plot.title = element_text(hjust = 0.5)) +
    guides(fill=FALSE) +
    labs(title=plotTitle) 
}




PlotSimplexAlternative=function(A,B,C,N){
  # can choose Uhat or U for Ufunction. If Uhat, need to Seed.
  # Seed(A,B)
  USimplex=UoverSimplex(A,B,C,N, U)

  Uround=signif(USimplex$U+.0000001,digits=4)
  USimplex$maximizer=(Uround==max(Uround))
  USimplex$fontf=sapply(USimplex$maximizer, function(l) (if (l) "bold" else "plain"))
  
  
  # cornerlables=data.frame(txt=sapply(1:3, function(i) paste("n", i, "=N", sep="")),
  #                         x=c(N+1, N/2, -1),
  #                         y=sqrt(.75)*c(-1, N+1, -1))
  
  cornerlables=data.frame(txt=sapply(1:3, function(i) paste("n", i, "=N", sep="")),
                          n1=c(N+1, -.5, -.5),
                          n2=c(-.5, N+1, -.5))
  


  plotTitle=bquote(alpha ~" = (" ~ 
                    .(paste(A, collapse=", ")) ~ 
                    "), "  ~ beta ~ " = ("  ~
                   .(paste(B, collapse=", ")) ~ 
                    ")")  
  ggplot(USimplex, aes(x=n1+sqrt(.25)*n2, y=sqrt(.75)*n2)) +
    geom_point(aes(color = U), size=sqrt(700/N))  + 
    scale_color_gradient( low = "skyblue4", high = "white") +
    geom_text(aes(label=format(U+.0000001, digits=3), #note addition of small number is to generate consistent rounding
                  fontface=fontf),
              size=10/N,
              show.legend=FALSE)  +
    scale_x_continuous(breaks = 0:N) + scale_y_continuous(breaks = 0:N) +
    geom_text(data=cornerlables,
              aes(label=txt),
              size= 15/N,
              show.legend=FALSE)  +
    theme_void() + 
    theme(legend.position="none", 
          panel.background = element_rect(fill = backcolor, colour = NA),
          plot.title = element_text(hjust = 0.5)) +
    coord_fixed(ratio = 1, xlim = c(-1.5,N+1.5), ylim=c(-1.5, sqrt(.75)*(N+1.5))) +
    labs(title=plotTitle)
}






SimplexPanel=function(N, alternativeplot=FALSE){
  C=c(0,0,0)
  AM=matrix(c(2,2,2,   #design matrix for plots
              2,2,3,
              2,2,1,
              3,3,2,
              3,3,1),
            5,3, byrow=TRUE)
  BM=matrix(4,5,3)-AM #corresponsing to uniform prior and pilot experiment with 2 units each
  
  for (i in 1:5) {
    if (!alternativeplot){
      PlotSimplex(AM[i,],BM[i,],C,N)
      filename=paste(c("../Figures/SimplexOrthogonal/ESWF_prior", AM[i,],"sample", N ,".pdf"), collapse="")
    } else {
      PlotSimplexAlternative(AM[i,],BM[i,],C,N)
      filename=paste(c("../Figures/SimplexTriangle/ESWF_prior", AM[i,],"sample", N ,"Alternative.pdf"), collapse="")          
    }
    ggsave(filename, width = 5, height = 5)
  }
  
}


OptimalPilot=function(A,B,C,M, parallel=TRUE){
  
  pilot=tibble(n1=0:M, Vn1=rep(0,M+1))
  
  if (parallel) {
    library(parallel)
    no_cores = detectCores()
    clust = makeCluster(no_cores, type="FORK")  #forking requires mac or Linux OS!
    #for Windows use this instead:
    #clust = makeCluster(no_cores)
    #clusterExport(clust, c("V", "UoverSimplex", "U", "SWF", "betabinomialvector", "betabinomial"))
    pilot$Vn1=parSapply(clust, pilot$n1, function(n)  V(A,B,C,c(n, M-n)))
    stopCluster(clust)
  } else {
    pilot$Vn1=sapply(pilot$n1, function(n)  V(A,B,C,c(n, M-n)))
  }
  
  ggplot(pilot, aes(x=as.integer(n1), y=Vn1)) +
    geom_line(size=1, color= "skyblue4")+
    xlab(expression(N[1])) + ylab(expression(V[0])) +
    scale_x_continuous(breaks = 0:M) +
    theme_light() +
    theme(panel.grid.minor = element_blank(),
          #panel.background = element_rect(fill = backcolor, colour = NA),
          axis.line.x = element_line(size = 0.5, colour = "black"))
}
