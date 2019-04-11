## Load libraries
library(shiny)
library(tidyverse)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  
  output$default <- renderText({ input$user_gene_clustering })
  
})
