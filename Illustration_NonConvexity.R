source("OptimalAssignmentFunctions/WelfareFunctions.R")
source("OptimalAssignmentFunctions/welfareplotsGraphics.R")
source("OptimalAssignmentFunctions/SimulatedWelfareFunctions.R")

source("IllustrationFunctions/Illustration_NonConvexity_Functions.R")

backcolor="azure2" #background color for plots
gridcolor="azure1"
fillcolor="skyblue4"

# expected welfare calculations

#2 treatment values
#de-facto degenerate prior for treatment 2
N=20
C=c(0,0)


A=c(1,10000)
B=c(5,10000)
stylizedDesign(A,B,C,N)

A=c(5,10000)
B=c(1,10000)
stylizedDesign(A,B,C,N)

A=c(3,10000)
B=c(3,10000)
stylizedDesign(A,B,C,N)



#and now for power calculations based on the same designs and prior means

powerCalc(c(1/6, .5),N)
powerCalc(c(5/6, .5),N)
powerCalc(c(.5, .5),N)


#and now for MSE of difference in means based on the same designs and prior means

MSEcalc(c(1/6, .5),N)
MSEcalc(c(5/6, .5),N)
MSEcalc(c(.5, .5),N)

