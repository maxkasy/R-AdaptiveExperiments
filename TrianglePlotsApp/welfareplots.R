rm(list = ls())

#setwd(dirname(parent.frame(2)$ofile))
source("welfareplotsFunctions.R")
source("welfareplotsGraphics.R")
source("SimulatedFunctions.R")

for (N in c(3,4,6,10)) {
  SimplexPanel(N, alternativeplot=TRUE)
}


k=3
A=rep(1,k)
B=rep(1,k)
C=rep(0,k)
M=10
OptimalPilot(A, B, C, M)

filename=paste(c("../../Figures/OptimalPilot/OptimalPilot_", M, "_prior", A ,".pdf"), collapse="")  
ggsave(filename, width = 5, height = 4)
  




