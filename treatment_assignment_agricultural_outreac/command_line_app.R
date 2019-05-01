##############################################
#Invoke this script from the terminal as follows:
#RScript --vanilla command_line_app.R priordata_test_nocovar.csv 50
#the first argument is the input file of prior data
#the second argument is the number of units in the current wave
##############################################

source("ReadDataApp.R")
source("modified_thompson.R")
source("optimal_stopping.R")

command_arguments = commandArgs(trailingOnly = T)
datapath = command_arguments[1]
Nt = as.integer(command_arguments[2])


priordata = ReadDataApp(datapath)


Dt = DtchoiceThompson_modified(priordata$Y,
                               priordata$D,
                               priordata$k,
                               Nt = Nt)


filename = paste(Sys.Date(), "_treatment_assignment.csv", sep = "")
write_csv(tibble(D = Dt), path = filename)

regret_reduction = expected_return(priordata$Y,
                                   priordata$D,
                                   priordata$k,
                                   Nt = Nt)

print("Expected reduction of regret: ")
print(regret_reduction)
