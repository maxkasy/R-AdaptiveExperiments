stylizedDesign=function(A,B,C, N){
  USimplex=UoverSimplex(A,B,C,N, U)
  #USimplex$deltaU=USimplex$U-c(USimplex$U[1],USimplex$U[1:N])
  
  ggplot(USimplex[-1,], aes(x=n1, y=U)) +
    geom_line(size=1, color= fillcolor)+
    xlab(expression(n[1])) + ylab("ESWF") +
    scale_x_continuous(breaks = seq(0,N,5)) +
    theme_light() +
    theme(panel.grid.minor = element_blank(),
          #panel.grid.major = element_line(colour=gridcolor),
          #panel.background = element_rect(fill = backcolor, colour = NA),
          axis.line.x = element_line(size = 0.5, colour = "black"))
    
  
  filename=paste(c("../Figures/ValueofInfo/ESWF_prior", A[1], B[1],"_stylized.pdf"), collapse="") 
  
  ggsave(filename, width = 5, height = 4)
}
 

#and now for power calculations based on the same designs and prior means
powerCalc=function(theta, N){
  
  #power of t-test based on normal approximation
  n1=1:N
  stdEst=sqrt(theta[1]*(1-theta[1])/n1)
  meanEst= theta[1]-theta[2]
  meanTstat=meanEst/stdEst
  critval=qt(.975,df=n1-1)
  power=pt(-critval-meanTstat, df=n1-1)+1-pt(critval-meanTstat, df=n1-1)
  power[1]=0 #no power for degrees of freedom less than 1
  power[2]=0
  
  powerData=data.frame(n1=n1, power=power)
  
  ggplot(powerData, aes(x=n1, y=power)) +
    geom_line(size=1, color= "skyblue4")+
    xlab(expression(n[1])) + ylab("power") +
    scale_x_continuous(breaks = seq(0,N,5)) +
    theme_light() +
    theme(panel.grid.minor = element_blank(),
          #panel.grid.major = element_line(colour=gridcolor),
          #panel.background = element_rect(fill = backcolor, colour = NA),
          axis.line.x = element_line(size = 0.5, colour = "black"))
  
  filename=paste(c("../Figures/ValueofInfo/powerCalc", round(100*theta[1]),"_stylized.pdf"), collapse="") 
  
  ggsave(filename, width = 5, height = 4)
}



#and now for MSE of difference in means based on the same designs and prior means
MSEcalc=function(theta, N){
  
  #power of t-test based on normal approximation
  n1=1:N
  MSE=theta[1]*(1-theta[1])/n1
  
  powerData=data.frame(n1=n1, MSE=MSE)
  
  ggplot(powerData, aes(x=n1, y=-MSE)) +
    geom_line(size=1, color= "skyblue4")+
    xlab(expression(n[1])) + ylab("-MSE") +
    scale_x_continuous(breaks = seq(0,N,5)) +
    theme_light() +
    theme(panel.grid.minor = element_blank(),
          #panel.grid.major = element_line(colour=gridcolor),
          #panel.background = element_rect(fill = backcolor, colour = NA),
          axis.line.x = element_line(size = 0.5, colour = "black"))
  
  filename=paste(c("../Figures/ValueofInfo/MSECalc", round(100*theta[1]),"_stylized.pdf"), collapse="") 
  
  ggsave(filename, width = 5, height = 4)
}

