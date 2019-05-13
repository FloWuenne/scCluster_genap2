#source("../uifunctions.R")
library(DT)
library(shinycssloaders)
library(shiny)



# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Visualize your scRNA-seq clustering results"),
  br(),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      textInput("user_gene_clustering", label = "Enter Genesymbol", value = "GAPDH"),
      
      ## Dropdown menu for color panel
      selectInput(
        inputId="color_palette", label="Choose a color palette for expression!", 
        choices=c("A","B","C","D","E"),
        selected="D", multiple=FALSE, selectize=FALSE),
      
      ## Scale menu for point size 
      sliderInput(
        inputId="point_size", label="Change cell point size!", 
        min =1, max = 6, value = 3, step = 1, round = TRUE) ,
      
      h4("Do you want to show cell type labels?"),
      checkboxInput("show_labels",
                    label = "Show labels",
                    value = TRUE)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("tSNE Clustering",
                 plotOutput("tsne_plot_cluster"),
                 plotOutput("tsne_plot_gene_expression")
        
      ),
      
      tabPanel("Marker",
               dataTableOutput("table_marker_genes")
      )
    
  )
  
  )
)
))
