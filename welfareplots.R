#setwd(dirname(parent.frame(2)$ofile))

library(tidyverse)

source("OptimalAssignmentFunctions/WelfareFunctions.R")
source("OptimalAssignmentFunctions/welfareplotsGraphics.R")
source("OptimalAssignmentFunctions/SimulatedWelfareFunctions.R")

Simplices = F #whether to plot welfare simplices
WaveDivision = T #whether to explore effect of wave sizes on welfare

if (Simplices) {
    for (N in c(3,4,6,10)) {
      SimplexPanel(N, alternativeplot=T)
    }
}


if (WaveDivision){
    k=3
    A=rep(1,k)
    B=rep(1,k)
    C=rep(0,k)
    M=10
    OptimalPilot(A, B, C, M)
    
    filename=paste(c("../Figures/OptimalPilot/OptimalPilot_", M, "_prior", A ,".pdf"), collapse="")  
    ggsave(filename, width = 5, height = 4)
}  




