#DirectoryList=list.dirs(recursive = F)
library(purrr)
library(rmarkdown)

listFunctions <- function(filename) {
  temp.env <- new.env()
  sys.source(filename, envir = temp.env)
  functions <- lsf.str(envir=temp.env)
  rm(temp.env)
  return(functions)
}

printfunctions=function(script) {
  cat("###", script, "\n")
  print(listFunctions(script))
  cat("\n\n\n")
}

printfunctions_insubfolders = function(){
  dirlist= dir(recursive=F, pattern="*Functions")
  for (directory in dirlist) {
    setwd(directory)
    cat("# ", directory[[1]], "\n\n")
    ScriptList=list.files(pattern="*.R", recursive = F)
    walk(ScriptList, printfunctions)
    setwd("../")
  }
}

sink("functionlist.Rmd")
printfunctions_insubfolders ()
sink()


rmarkdown::render("functionlist.Rmd", "pdf_document")









  
