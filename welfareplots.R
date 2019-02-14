#setwd(dirname(parent.frame(2)$ofile))

library(tidyverse)

source("OptimalAssignmentFunctions/WelfareFunctions.R")
source("OptimalAssignmentFunctions/welfareplotsGraphics.R")
source("OptimalAssignmentFunctions/TreatmentAssignmentFunctions.R")
source("OptimalAssignmentFunctions/SimulatedWelfareFunctions.R")

Simplices = F #whether to plot welfare simplices
OptimalAssignments = F #whether to plot optimal assignments
WaveDivision = F #whether to explore effect of wave sizes on welfare
ThompsonMapping = T

if (Simplices) {
    for (N in c(3,4,6,10)) {
      SimplexPanel(N, alternativeplot=T)
    }
}

if (OptimalAssignments) {
    for (N2 in c(3,4,6,10)) {
        PlotOptimalAssignment(n1=rep(2,3), N2)
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


if (ThompsonMapping) {
    p_list=list(a=c(.2,.2333, .2666, .3),
                b=c(.1,.2,.3,.4),
                c=c(.05,.1, .4,.45),
                d=c(.05,.1,.15,.7))
    ThompsonMappingPlots(p_list)
}

