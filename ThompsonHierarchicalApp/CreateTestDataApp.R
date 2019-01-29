library(tidyverse)

rr=1000

CreateTestData=function() {
    priordata=tibble(
        covar1=rep(1:4, 3*rr),
        treatment1=rep(c(1,0,0), 4*rr),
        treatment2=rep(c(0,1,0), 4*rr),
        treatment3=rep(c(0,0,1), 4*rr),
        theta_i=.1* covar1 + .1 * treatment1  + .2 * treatment2 + .6*treatment3 - .2*covar1*treatment3,
        outcome=runif(12*rr)<theta_i
    )
    
    
    newwave=tibble(
      covar1=rep(1:4, 3*rr)
    )
    
    write_csv(priordata, "priordata_test_v2.csv")
    write_csv(newwave, "newwave_test_v2.csv")
}

CreateTestData()