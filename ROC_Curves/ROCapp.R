#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(MutationalPatterns)
library(ROCR)
source("bayes_function.R")
source("simulation.R")
source("get_classification_df.R")
source("get_performance_object.R")

ffpe.signature <- as.matrix(load_old_mutational_matrix("supplied_results/ffpe.signature.txt"))

#create a vector of the 96 snv mutation types:
##this is used to assign rownames to COSMIC signature matrices later on
mutations <- rownames(ffpe.signature)

#load COSMIC version 3 signatures
cosmic.v3 <- get_known_signatures(muttype = "snv")

cosmic <- colnames(cosmic.v3)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Bayesian Classification of COSMIC and FFPE mutations"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput(inputId = "total_mutations",
                     label = "Total number of mutations:",
                     min = 100,
                     max = 1000,
                     value = 100),
         sliderInput(inputId = "proportion_ffpe",
                     label = "Proportion of total mutations from the FFPE signature",
                     min = 0,
                     max = 1,
                     value = 0.5),
         selectInput(inputId = "signature",
                      label = "COSMIC Signature:",
                      choices = cosmic,
                      selected = "SBS1")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("plot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input,output) {
  
  
  performance.object <- reactive(get_performance_object(input$signature, input$total_mutations,
                                                        input$proportion_ffpe))
  
  
  output$plot <- renderPlot({
    p <- plot(performance.object(), avg = "threshold")
    p
  })

}

# Run the application 
shinyApp(ui = ui, server = server)

