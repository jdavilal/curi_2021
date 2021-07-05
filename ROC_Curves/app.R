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

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Bayesian Classification of COSMIC and FFPE mutations"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("total_mutations",
                     "Total number of mutations:",
                     min = 100,
                     max = 1000,
                     value = 100),
         sliderInput("proportion_ffpe",
                     "Proportion of total mutations from the FFPE signature",
                     min = 0,
                     max = 1,
                     value = 0.5),
         selectInput(inputId = "signature",
                     label = "COSMIC Signature:",
                     choices = c ("SBS1", "SBS2", "SBS3","SBS4", "SBS5", "SBS6",
                                  "SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS8",
                                  "SBS9", "SBS10a", "SBS10b", "SBS10c", "SBS10d",
                                  "SBS11", "SBS12", "SBS13", "SBS14", "SBS15",
                                  "SBS16", "SBS17a", "SBS17b", "SBS18", "SBS19",
                                  "SBS20", "SBS21", "SBS22", "SBS23", "SBS24", "SBS25",
                                  "SBS26", "SBS28", "SBS29", "SBS30", "SBS31",
                                  "SBS32", "SBS33", "SBS34", "SBS35", "SBS36", "SBS37",
                                  "SBS38", "SBS39", "SBS40", "SBS41", "SBS42", "SBS44",
                                  "SBS84", "SBS85", "SBS86", "SBS87", "SBS88", "SBS89",
                                  "SBS90", "SBS91", "SBS92", "SBS93", "SBS94"),
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
  
  ffpe.signature <- as.matrix(load_old_mutational_matrix("supplied_results/ffpe.signature.txt"))
  
  #create a vector of the 96 snv mutation types:
  ##this is used to assign rownames to COSMIC signature matrices later on
  mutations <- rownames(ffpe.signature)
  
  #load COSMIC version 3 signatures
  cosmic.v3 <- get_known_signatures()
  
  #extract signature of interest from the COSMIC v3 matrix
  #cosmic.signature <- cosmic.v3[,str(input$signature)]
  #assign rownames to the signature matrix
  #rownames(cosmic.signature) <- mutations
  
  #cosmic_sig_name <- reactive({
    #req(input$signature)
    #cosmic_sig_name <- input$signature
  #})
  cosmic_sig_name <- reactive(input$signature)
  
  cosmic.signature1 <- cosmic.v3[,cosmic_sig_name]
  cosmic.signature <- as.matrix(cosmic.signature1)
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
  #cosmic.sample <- create_signature_sample_vector(cosmic.signature, 50)
  ffpe.sample <- create_signature_sample_vector(ffpe.signature, ffpe.mutation.number)
  #ffpe.sample <- create_signature_sample_vector(ffpe.signature, 50)
  
  
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

