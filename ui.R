#source("../uifunctions.R")
library(DT)
library(shinycssloaders)
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(plotly)
library(shinyFiles)
library(shinydashboard)

## Define how much RAM R should ask for for caching
options(shiny.maxRequestSize = 10000*1024^2)

# Define UI for application that draws a histogram
shinyUI(
  
  fluidPage(theme = shinytheme("readable"),
            
            tags$head(tags$style(
              type="text/css",
              "#image img {max-width: 100%; width: auto; height: auto}"
            )),
            
            
            ## include css stylesheet
            #includeCSS("styles.css"),
            
            id = "main_page",
            
            # Application title
            h1("Visualize and analyse your scRNA-seq clustering results",
                       windowTitle = "scRNA-seq clustering"),
            
            hr(),
            
            fluidRow(style = "background-color:#C2C5CC;",align="center",
                     br(),
                     column(4,align="center",
                            conditionalPanel(condition = "output.user_cluster_labels",
                                      uiOutput("available_cluster_labels"))
                            ),
                     column(4,align="center",
                            conditionalPanel(condition = "output.user_cluster_labels",
                                             icon("book","fa-3x",  lib = "font-awesome"))
                     ),
                     column(4,align="center",
                            conditionalPanel(condition = "output.user_cluster_labels",
                                             p("You are currently using the following annotation:")
                                             ),
                            conditionalPanel(condition = "output.user_cluster_labels",
                                             h4(textOutput("selected_annotation"))
                                             )
                            ),

                     br()
                     ),
            
            hr(),
  
  tabsetPanel(
    
    tabPanel("Data upload", icon = icon("file-upload"),

             br(),
             br(),
             
             fluidRow(
               column(12, align="center",
                      h1("Upload files"),
                      
                      shinyDirButton("file_dir", "Folder select", "Please select a folder"),
                      
                      verbatimTextOutput("test_path")
                      )
               ), 
             
             fluidRow(
               
               hr(),
               
               column(12, align="center",
                      imageOutput("genap_logo")
                      #img(src='GenAP_powered_reg.png', align = "center")
               )
             )


    ),
    
    ## Table with dimensional reduction
    tabPanel("Dimensional reduction",icon = icon("circle"),
             
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
             conditionalPanel(condition = "!output.dimredoutput",
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
    tabPanel("Marker genes",icon = icon("highlighter"),
             
             ## Show this text only when no marker table has been uploaded!
             conditionalPanel(condition = "!output.marker_genes_table_out",
                              h3("Please upload a marker gene table!")
             ),
             
             conditionalPanel(condition = "output.marker_genes_table_out",
                              dataTableOutput("table_marker_genes")
             )
            ),
    
    ## Panel for renaming clusters
    tabPanel("Modify & add annotations",icon = icon("file-signature"),
             
             br(),

             fluidRow(
               
               conditionalPanel(condition = "!output.dimredoutput",
                                h3("Please upload a folder with processed files with clustering results!")
               ),
               
               column(4,
                      conditionalPanel(condition = "output.dimredoutput",
                      textInput(inputId = "user_added_cluster", 
                                label = "Please enter a name for your new annotation?",
                                value = "")
                      )),
               
               column(4,
                      conditionalPanel(condition = "output.dimredoutput",
                                       actionButton("add_annotation", "Add new annotation!"),
                                       br()
                      ))
             ),
             
             br(),
             hr(),
             
             
             fluidRow(
               column(4,

                      conditionalPanel(condition = "output.dimredoutput",
                                       selectInput("rename_method", 
                                                   label = "How do you want to rename your clusters:", 
                                                   choices = c("Assigned clusters" = "assigned_clusters",
                                                               "Gene expression" = "gene_expression",
                                                               "Cell selection" = "cell_selection"),
                                                   selected = "Assigned clusters", 
                                                   multiple = FALSE)
                                       )
                      ),
               
               column(6,
                      conditionalPanel(condition = "output.rename_selected == 'assigned_clusters' && output.dimredoutput",
                                       p("This panel let's you rename clusters based on previously assigned groups of cell 
                                         clusters!")),
                      
                      conditionalPanel(condition = "output.rename_selected == 'gene_expression' && output.dimredoutput",
                                       p("This panel let's you rename cell clusters based on gene expression thresholds.
                                         This can be a powerful way of classifying cells into clusters of shared expression
                                         profiles based on known markers. Currently supports up to 2 markers!"))
                      )
               
             ),
             
             br(),
             hr(),
             
             
             ## Reactive panels that change based on how the user wants to rename clusters
             fluidRow(
               column(4,
                      conditionalPanel("output.rename_selected == 'assigned_clusters'",
                                       plotOutput("dimred_plot_rename_assigned_clusters")
                                       ),
                      
                      conditionalPanel("output.rename_selected == 'gene_expression'",
                                       plotOutput("dimred_plot_rename_expression")
                      ),
                      
                      conditionalPanel("output.rename_selected == 'gene_expression'",
                                       plotOutput("rename_expression_dens", height = "200px")
                      )

                      ),
               
               column(4,
                      uiOutput("rename_list"),
                      
                      uiOutput("gene_exp_threshold"),
                    
                      uiOutput("plot_gene_rename_button"),
                      br(),
                      br(),
                      h4(textOutput("cells_exp_selected")),
                      
                      tableOutput("test")
                      ),
               
               column(4,
                      uiOutput("cluster_anno_text"),
                      br(),
                      uiOutput("save_anno_button")
                      )
      
            
             # 
             # verbatimTextOutput("brush")

             )
    )
                     

    
  ) # end of TabsetPanel
  ) # end of fluidPage
) # End of shinyUI
