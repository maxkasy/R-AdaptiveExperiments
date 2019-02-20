source("OptimalAssignmentFunctions/WelfareFunctions.R")
source("OptimalAssignmentFunctions/TreatmentAssignmentFunctions.R")
source("OptimalAssignmentFunctions/welfareplotsGraphics.R")
source("OptimalAssignmentFunctions/SimulatedWelfareFunctions.R")

source("CalibratedSimulationFunctions/ReadData.R")
source("CalibratedSimulationFunctions/calibratedSimulationFunctions.R")

source("ThompsonHierarchicalFunctions/MCMC_HierarchicalThompson.R")
source("ThompsonHierarchicalFunctions/calibratedSimulationFunctionsCovariates.R")

#backcolor="azure2" #background color for plots
#gridcolor="azure1"
fillcolor="skyblue4"
today=Sys.Date()

printFigures=T
DoNoCovariates=T
DoHistograms=T
DoCovariates=T

DataList=ReadAllData(printFigures)

#number of replications
MC_replicates=20000


#Simulations without covariates
if (DoNoCovariates) {
    #which methods of treatment assignment to simulate
    methods=c(12,9,8,1)
    #multiples of original sample size
    multiplicities=c(.5, 1, 1.5)
    #number of waves
    wavenumbers=c(2,4,10)
    #columnames for tables
    columnames=paste(wavenumbers, "waves")
    
    for (Mult in multiplicities) {
      for (application in 1:length(DataList)) {
        applicationData=rep(DataList[application], length(wavenumbers))
        for (i in 1:length(wavenumbers)) {
          applicationData[[i]]$wavesizes = rep(floor(applicationData[[i]]$N * Mult / wavenumbers[i]) ,wavenumbers[i])
          applicationData[[i]]$waves=wavenumbers[i]
        }
        filename=paste(today,"_CalibratedSimulations_", DataList[[application]]$filename, "_" ,Mult, sep="")
        DesignTable(applicationData,methods,MC_replicates,columnames,filename, DoHistograms)
      }
    }
}


MC_replicates=5000

#Simulations with covariates
if (DoCovariates) {
    #which methods of treatment assignment to simulate
    methods=c(5,4,3,2,1)
    #multiples of original sample size
    multiplicities=c(.5, 1, 1.5)
    # number of waves
    wavenumbers=c(2,4,10)
    #columnames for tables
    columnames=paste(wavenumbers, "waves")
    
    for (Mult in multiplicities) {
      for (application in 1:length(DataList)) {
        applicationData=rep(DataList[application], length(wavenumbers))
        for (i in 1:length(wavenumbers)) {
          applicationData[[i]]$wavesizes = rep(floor(applicationData[[i]]$N * Mult / wavenumbers[i]) ,wavenumbers[i])
          applicationData[[i]]$waves=wavenumbers[i]
        }
        filename=paste(today,"_Covariates_CalibratedSimulations_", DataList[[application]]$filename, "_" ,Mult, sep="")
        DesignTableCovariates(applicationData,methods,MC_replicates,columnames,filename)
      }
    }
}









