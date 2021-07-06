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

ffpe.signature <- as.matrix(load_old_mutational_matrix("supplied_results/ffpe.signature.txt"))

#create a vector of the 96 snv mutation types:
##this is used to assign rownames to COSMIC signature matrices later on
mutations <- rownames(ffpe.signature)

#load COSMIC version 3 signatures
cosmic.v3 <- get_known_signatures()

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
                      choices = c ("SBS1" = 1, "SBS2" = 2, "SBS3" = 3,"SBS4" = 4, "SBS5" = 5, "SBS6" = 6,
                                   "SBS7a" = 7, "SBS7b" = 8, "SBS7c" = 9, "SBS7d" = 10, "SBS8" = 11,
                                   "SBS9" = 12, "SBS10a" = 13, "SBS10b" = 14, "SBS10c" = 15, "SBS10d" = 16,
                                   "SBS11" = 17, "SBS12" = 18, "SBS13" = 19, "SBS14" = 20, "SBS15" = 21,
                                   "SBS16" = 22, "SBS17a" = 23, "SBS17b" = 24, "SBS18" = 25, "SBS19" = 26,
                                   "SBS20" = 27, "SBS21" = 28, "SBS22" = 29, "SBS23"= 30 , "SBS24" =31, 
                                   "SBS25" = 32,"SBS26" = 33, "SBS28" = 34, "SBS29" = 35, "SBS30" = 36, 
                                   "SBS31" = 37, "SBS32" = 38, "SBS33" = 39, "SBS34" = 40, "SBS35" = 41, 
                                   "SBS36" = 42, "SBS37" = 43, "SBS38" = 44, "SBS39" = 45, "SBS40" = 46, 
                                   "SBS41" = 47, "SBS42" = 48, "SBS44" = 49, "SBS84" = 50, "SBS85" = 51,
                                   "SBS86" = 52, "SBS87" = 53, "SBS88" = 54, "SBS89" = 55,
                                   "SBS90" = 56, "SBS91" = 57, "SBS92" = 58, "SBS93" = 59, "SBS94" = 60),
                      selected = "SBS1")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("plot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  #extract signature of interest from the COSMIC v3 matrix
  #assign rownames to the signature matrix
  
  cosmic_sig_name <- reactive({
    req(input$signature)
    cosmic_sig_name <- input$signature
  })

  cosmic.signature <- as.matrix(cosmic.v3[,cosmic_sig_name])
  rownames(cosmic.signature) <- mutations
  
  ffpe.mutation.number <- reactive({
    req(input$total_mutations)
    req(input$proportion_ffpe)
    ffpe.mutation.number <- as.double(input$total_mutations*input$proportion_ffpe)
  })
  
  cosmic.mutation.number <- reactive({
    req(input$total_mutations)
    req(input$proportion_ffpe)
    cosmic.mutation.number <- as.double(input$total_mutations-(input$total_mutations*input$proportion_ffpe))
  })
  
  #create sample vectors for COSMIC and FFPE mutations
  cosmic.sample <- create_signature_sample_vector(cosmic.signature, cosmic.mutation.number)
  ffpe.sample <- create_signature_sample_vector(ffpe.signature, ffpe.mutation.number)
  
  
  class.df <- get_classification_df(list(cosmic.sample, ffpe.sample), c(cosmic_sig_name, "FFPE"),
                                    list(cosmic.signature, ffpe.signature))
  
  class.df2 <- class.df %>%
    mutate(cosmic.indicator = ifelse(truth == cosmic_sig_name, 1, 0))
  
  
  predictions <- class.df2$cosmic_sig_name
  labels <- class.df2$cosmic.indicator
  pred <- prediction(predictions, labels)
  
  perf <- performance(pred, "tpr", "fpr")
  
  output$plot <- renderPlot({
    
     p <- plot(perf,
          avg="threshold")
     p
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

