# to upload, run rsconnect::deployApp('ThompsonHierarchical')

library(shiny)
rm(list = ls())
source("ThompsonHierarchical.R")
source("ReadDataApp.R")


ui <- fluidPage(
#  titlePanel(h1("Optimal treatment assignment given covariates", align = "center")),
  verticalLayout(
    includeMarkdown("instructions.md"),
    hr(),
    fluidRow(column(4, fileInput("file1", "Choose CSV file of previous data",
                              multiple = FALSE,
                              accept = c("text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv"))),
             column(4, fileInput("file2", "Choose CSV file of new wave",
                              multiple = FALSE,
                              accept = c("text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv"))),
             column(4,numericInput(inputId="RR", label="Replicate draws", value=10, min=1, max=100))
             ),
    fluidRow(column(4,actionButton(inputId = "calcbutton", label = "Calculate treatment assignment")),
             column(4,downloadButton("downloadData", "Download treatment assignment"))),
    hr(),
    
    fluidRow(column(6,tableOutput("treatmentcounts")),
             column(6,tableOutput("designtable")))
  )
)



server <- function(input, output, session) {
  
  v = reactiveValues()
  
  observeEvent(input$calcbutton,{
    req(input$file1)
    req(input$file2)
    
    #loading files
    #make sure that options are chosen right!
    priordata=ReadDataApp(input$file1$datapath)    
    newwave=ReadDataApp(input$file2$datapath,
                        priordata$key)    
    
    #calculating treatment assignment
    if (priordata$nx > 1) {
      newwave$Dstar=as.integer(
                  DtchoiceThompsonHierarchicalModified(priordata$Y,priordata$D,priordata$X, #outcomes, treatments, and covariates thus far
                                                priordata$k,priordata$nx, #number of treatments and number of strata
                                                newwave$Xt, RR=input$RR))
      v$newwave =rename(newwave,
                        stratum =Xt,
                        treatment=Dstar)
      
      v$treatmentcounts=    v$newwave %>%
        group_by(stratum, treatment) %>%
        summarise(count=n()) %>%
        spread(treatment, count)
    } else {
      newwave$Dstar=as.integer(
                    DtchoiceThompsonModified(priordata$Y,priordata$D, #outcomes and treatments thus far
                                     priordata$k, #number of treatments
                                     nrow(newwave),
                                     RR=input$RR))
      v$newwave =rename(newwave,
                        treatment=Dstar)
      
      v$treatmentcounts=    v$newwave %>%
        group_by(treatment) %>%
        summarise(count=n())
    }
    

    

  })
  
  output$designtable =  renderTable(
    {v$newwave}, 
    align="c", 
    caption="Treatment assignment",
    caption.placement = getOption("xtable.caption.placement", "top")
   )
  
  output$treatmentcounts =  renderTable(
    {v$treatmentcounts}, 
    align="c",
    caption="Treatment counts",
    caption.placement = getOption("xtable.caption.placement", "top")
  )
 
#download optimal design
  output$downloadData <- downloadHandler(
    filename = "treatmentassignment.csv",
    content = function(file) {
      write_csv(v$newwave, file)
    }
  )
  

 
}

# Run the app ----
shinyApp(ui = ui, server = server)