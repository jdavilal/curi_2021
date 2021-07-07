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
library(gridExtra)
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
         tabsetPanel(type = "tabs",
                     tabPanel("ROC Curve", plotOutput("rocplot"), textOutput("auc")),
                     tabPanel("Reconstructed Profile", plotOutput("reconplot"), textOutput("cos_sim"))
         )
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input,output) {
  
  
  performance.object <- reactive(get_performance_object(input$signature, input$total_mutations,
                                                        input$proportion_ffpe))
  
  profile.matrices <- reactive(get_profile_matrices(input$signature, input$total_mutations,
                                           input$proportion_ffpe))
  
  output$rocplot <- renderPlot({
    p <- plot(performance.object()[[1]], avg = "threshold")
    p
  })
  
  output$auc <- renderText({
    paste0("Area Under Curve: ", performance.object()[[2]])
  })
  
  
  output$reconplot <- renderPlot({
    p1 <- plot_96_profile(profile.matrices()[[1]])+
      labs(title = "True COSMIC Signature Profile")
    #plot_96_profile for COSMIC sample vector
    p2 <- plot_96_profile(profile.matrices()[[2]])+
      labs(title = "Mutational Profile of COSMIC Sample")
    #plot_96_profile for COSMIC and FFPE sample vectors
    p3 <- plot_96_profile(profile.matrices()[[3]])+
      labs(title = "Mutational Profile of COSMIC Sample with FFPE Added")
    #plot_96_profile for reconstructed COSMIC profile
    p4 <- plot_96_profile(profile.matrices()[[4]])+
      labs(title = "Reconstructed Profile After Removing FFPE According to Bayesian Classifier")
    grid.arrange(p1,p2,p3, p4)
  })
  
  output$cos_sim <- renderText({
    paste0("Cosine Similarity Between Original COSMIC Sample Profile and Reconstructed Profile: ", profile.matrices()[[5]])
  })

}

# Run the application 
shinyApp(ui = ui, server = server)

