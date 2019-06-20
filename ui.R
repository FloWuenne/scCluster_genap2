#source("../uifunctions.R")
library(DT)
library(shinycssloaders)
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(plotly)
library(shinyFiles)

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
    tabPanel("Marker genes",
             
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
             
             br(),
             fluidRow(
               column(4,
                      ## Show this text only when no marker table has been uploaded!
                      conditionalPanel(condition = "input.feather_file == 0",
                                       h3("Please upload a feather file with clustering results!")
                      ),
                      
                      selectInput("rename_method", 
                                  label = "How do you want to rename your clusters:", 
                                  choices = c("Assigned clusters" = "assigned_clusters",
                                              "Gene expression" = "gene_expression",
                                              "Cell selection" = "cell_selection"),
                                  selected = "Assigned clusters", 
                                  multiple = FALSE)),
               column(6,
                      conditionalPanel(condition = "output.rename_selected == 'assigned_clusters'",
                                       p("This panel let's you rename clusters based on previously assigned clusters by the clustering
                                         algorithm you used. This works best if you think the cluster assigned are valid and you 
                                         want to give them meaningful names!")),
                      
                      conditionalPanel(condition = "output.rename_selected == 'gene_expression'",
                                       p("This panel let's you rename cell clusters based on gene expression thresholds.
                                         This can be a powerful way of classifying cells into clusters of shared expression
                                         profiles based on known markers and expertise. Currently only supports 1 threshold!"))
                      
                      )
                      
             ),

             br(),
             hr(),
             
             fluidRow(
               column(4,
                      conditionalPanel("output.rename_selected == 'assigned_clusters'",
                                       plotOutput("dimred_plot_rename_assigned_clusters")
                                       ),
                      
                      conditionalPanel("output.rename_selected == 'gene_expression'",
                                       plotOutput("dimred_plot_rename_expression")
                      ),
                      
                      conditionalPanel("output.rename_selected == 'gene_expression'",
                                       plotOutput("rename_expression_hist", height = "200px")
                      )

                      ),
               
               column(4,
                      uiOutput("rename_list"),
                      
                      uiOutput("gene_exp_threshold"),
                      
                      uiOutput("plot_gene_rename_button"),
                      br(),
                      br(),
                      h4(textOutput("cells_exp_selected"))
                      ),
               column(4,
                      p("put the controls for renaming clusters here!"))
               
             )
             
             # 
             # verbatimTextOutput("brush")

    )
                     

    
  ) # end of TabsetPanel
  ) # end of fluidPage
) # End of shinyUI
