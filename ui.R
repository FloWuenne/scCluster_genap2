#source("../uifunctions.R")
library(DT)
library(shinycssloaders)
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(plotly)


# Define UI for application that draws a histogram
shinyUI(
  fluidPage(theme = shinytheme("cosmo"),
  
  # Application title
  titlePanel("Visualize and analyse your scRNA-seq clustering results",
             windowTitle = "scRNA-seq clustering"),
  br(),
  
  tabsetPanel(
    
    tabPanel("Dimensional reduction",
             
             br(),
             
             fluidRow(
               
               ## Dropdown menu for color panel
               column(width= 4,
                      selectInput(inputId="dimred_color_palette", label="Choose a color palette for expression!",
                                  choices=c("A","B","C","D","E"),
                                  selected="D", multiple=FALSE, selectize=FALSE) ),
               
               ## Scale menu for point size 
               column(width= 4,
                      sliderInput(inputId="point_size", label="Change cell point size!",
                                  min = 2, max = 10, value = 6, step = 2, round = TRUE) ),
               
               column(width= 4,
                      textInput("user_gene_clustering", label = "Enter Genesymbol", "GAPDH"),
                      actionButton("plot_gene_button", "Plot gene!") )
               ),
                 
                 br(),
                 br(),
             
             fluidRow(
               column(width=6,
                      plotlyOutput("dimred_plot_plotly") %>%
                        withSpinner(type = 8)
                      ),
               
               column(width=6,
                      plotlyOutput("dimred_gene_plot_plotly") %>%
                        withSpinner(type = 8)
                      )
               ),
             
             ## Violin plot and table
             fluidRow(
               
               ## Violin plot for gene selected
               column(width=12,
                      plotlyOutput("vlnplot_user_gene_plotly") %>%
                        withSpinner(type = 8)
                      )
               )
             
             ## old plots based on ggplot
             #withSpinner(plotOutput("tsne_plot_cluster")),
             #plotOutput("tsne_plot_gene_expression"),
             #plotOutput("vlnplot_user_gene")
        
      ),
    
    tabPanel("Marker_genes",
             dataTableOutput("table_marker_genes")
            )
    
  ) # end of TabsetPanel
  ) # end of fluidPage
) # End of shinyUI
