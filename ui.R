#source("../uifunctions.R")
library(DT)
library(shinycssloaders)
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(plotly)

options(shiny.maxRequestSize = 10000*1024^2)


# Define UI for application that draws a histogram
shinyUI(
  fluidPage(theme = shinytheme("cosmo"),
            id = "main_page",
            
            # Application title
            titlePanel("Visualize and analyse your scRNA-seq clustering results",
                       windowTitle = "scRNA-seq clustering"),
            
            br(),
  
  tabsetPanel(
    
    tabPanel("Data upload",
             

             br(),
             
             # # Ask user to upload their clustering results in feather format
             # fileInput("user_cluster_rds", "Choose a feather file containing clustering results...",
             #           multiple = FALSE),
             
             h3("Clustering feather file"),
             shinyFilesButton(id = "feather_file", 
                              label = "File select", 
                              title = "Please select a file", multiple = FALSE),
             
             #verbatimTextOutput("filepaths"),
             br(),
             
             h3("Gene names file"),
             shinyFilesButton(id = "gene_names", 
                              label = "File select", 
                              title = "Please select a file", multiple = FALSE),
             
             #verbatimTextOutput("filepaths"),
             
             h3("Marker genes table"),
             
             shinyFilesButton(id = "marker_genes", 
                              label = "File select", 
                              title = "Please select a file", multiple = FALSE)
    ),
    
    ## Table with dimensional reduction
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
             
             ## Show this text only when no feather file has been uploaded
             conditionalPanel(condition = "input.feather_file == 0",
                              h3("Please upload a feather file with clustering results!")
             ),
             
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
    
    ## Table with marker genes
    tabPanel("Marker_genes",
             
             ## Show this text only when no marker table has been uploaded!
             conditionalPanel(condition = "input.marker_genes == 0",
                              h3("Please upload a marker gene table!")
             ),
             
             conditionalPanel(condition = "input.marker_genes",
                              dataTableOutput("table_marker_genes")
             )
            ),
    
    ## Panel for renaming clusters
    tabPanel("Assign cluster names",
                     
                     h4("This panel will contain options to rename clusters. Currently I have 3 methods planned.
                     1) Rename clusters based on assigned communities.
                     2) Assign clusters based on expression.
                        3) Assign clusters based on manually selecting cells."))
    
  ) # end of TabsetPanel
  ) # end of fluidPage
) # End of shinyUI
