# to upload, run rsconnect::deployApp('webApp1Wave')
# app available at https://maxkasy.shinyapps.io/webapp1wave/

library(shiny)
rm(list = ls())
source("welfareplotsFunctions.R")
source("welfareplotsGraphics.R")
source("SimulatedFunctions.R")



ui <- fluidPage(
  titlePanel(h1("Optimal treatment assignment in the second wave", align = "center")),
  sidebarLayout(
    sidebarPanel(
      numericInput("N", 
                   label=h3("Total units in wave 2"), 
                   value = 4,max=20, min=0),
        h3("Outcomes of wave 1"),
        h5("treatment 1", align = "center"),
        fluidRow(
          column(6, numericInput("n1", 
                              label="Units", 
                              value = 2,min=0)),
          column(6, numericInput("s1", 
                              label="Successes", 
                              value = 1,min=0))),
        h5("treatment 2", align = "center"),
        fluidRow(
          column(6, numericInput("n2", 
                                 label="Units", 
                                 value = 2,min=0)),
          column(6, numericInput("s2", 
                                 label="Successes", 
                                 value = 1,min=0))),
        h5("treatment 3", align = "center"),
        fluidRow(
          column(6, numericInput("n3", 
                                 label="Units", 
                                 value = 2,min=0)),
          column(6, numericInput("s3", 
                                 label="Successes", 
                                 value = 1,min=0)))
    ),
    mainPanel(
      plotOutput("grid", width = "100%")
      # textOutput("benchmark")
      # textOutput("twowavebench")
    )
  )  
)



server <- function(input, output, session) {

  #make sure there can be no more successes than observations in each treatment arm
  observe({
    updateNumericInput(session, "s1", max=input$n1)
  })
  observe({
    updateNumericInput(session, "s2", max=input$n2)
  })
  observe({
    updateNumericInput(session, "s3", max=input$n3)
  })
  
  ## generate simplex plot
  output$grid=renderPlot({
      A=c(input$s1,input$s2,input$s3) + rep(1,3)
      B=c(input$n1,input$n2,input$n3) - c(input$s1,input$s2,input$s3) + rep(1,3)
      C=rep(0,3)
      PlotSimplexAlternative(A,B,C,input$N)
  })
  
  # output$benchmark = renderText({ 
  #   NN=input$N + input$n1 + input$n2 + input$n3
  #   paste("Expected welfare of one-wave experiment with ", 
  #         NN,
  #         "observations: ",
  #         V(rep(1,3), rep(1,3), rep(0,3), NN))
  # })
  
  # output$twowavebench = renderText({ 
  #   N2=input$N
  #   N1= input$n1 + input$n2 + input$n3
  #   paste("Expected welfare of two-wave experiment with ", 
  #         N1, " and ", N2,
  #         "observations: ",
  #         V(rep(1,3), rep(1,3), rep(0,3), c(N1,N2)))
  # })  
}

# Run the app ----
shinyApp(ui = ui, server = server)