#source("../uifunctions.R")
library(DT)
library(shinycssloaders)
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(plotly)
library(shinyFiles)
library(shinydashboard)
library(shinyalert)

## Define how much RAM R should ask for for caching
options(shiny.maxRequestSize = 10000*1024^2)

# Define UI for application that draws a histogram
shinyUI(
  
  fluidPage(theme = shinytheme("readable"),
            
            useShinyalert(),  # Set up shinyalert
            
            tags$head(tags$style(
              type="text/css",
              "#image img {max-width: 100%; width: auto; height: auto}"
            )),
            
            
            ## include css stylesheet
            #includeCSS("styles.css"),
            
            id = "main_page",
            
            # Application title
            h1("Visualize and annotate your scRNA-seq clustering results",
                       windowTitle = "scRNA-seq clustering"),
            
            hr(),
            
            fluidRow(style = "background-color:#B0DFE5;",align="center",
                     br(),
                     column(1,align="center",
                            conditionalPanel(condition = "output.user_cluster_labels",
                                             icon("book","fa-3x",  lib = "font-awesome"))
                     ),
                    column(4,align="center",
                           conditionalPanel(condition = "output.user_cluster_labels",
                                            p("Current dataset:")
                                            ),
                           conditionalPanel(condition = "output.user_cluster_labels",
                                            h4(textOutput("test_file_dir"))
                           )
                           ) ,
                     column(4,align="center",
                            conditionalPanel(condition = "output.user_cluster_labels",
                                      uiOutput("available_cluster_labels"))
                            ),
                     column(3,align="center",
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
                      h1("Upload dataset"),
                      
                      shinyDirButton("file_dir", "Choose a dataset!", "Please select a dataset!")
                    
                      )
               ), 
             
             fluidRow(
               
               hr(),
               
               column(12, align="center",
                      imageOutput("genap_logo", height = "80%")
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
                      textInput(inputId = "user_gene_clustering", label = "Enter Genesymbol"),
                      actionButton("plot_gene_button", "Plot gene!") )
               ),
                 
                 br(),
                 br(),
             
             ## Show this text only when no feather file has been uploaded
             conditionalPanel(condition = "!output.dimredoutput",
                              h3("Please upload a processed dataset!")
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
    
    ## Panel for renaming clusters
    tabPanel("Modify & add annotations",icon = icon("file-signature"),
             
             br(),
             
             fluidRow(
               column(6,align = "left",
                      p("This panel allows you to add annotations to you clustering. You can add a new
                        annotation column, you can delete existing columns and you can save annotations
                        for future use. When you create a new annotation column, the labels of the currently
                        selected annotation will initally be copied.")
                      ),
               column(6,align = "left",
                      p("You can group cells into clusters
                        using a combination of 3 different methods:
                        1) rename existing clusters
                        2) assign cells based on expression threshold of candidate genes
                        3) manually select a group of cells that you want to assign a cluster")
                      )
             ),

             fluidRow(
               column(12,
                      conditionalPanel(condition = "!output.dimredoutput",
                                       h3("Please upload a processed dataset!")
                                       )
                      )
             ),
             br(),
             
             fluidRow(
               
               column(3,align = "left",
                      conditionalPanel(condition = "output.dimredoutput",
                      h4(textInput(inputId = "user_added_cluster", width = "100%",
                                label = "Enter a name for your new annotation:",
                                value = ""))
                      )
                      ),
               
               column(3,align = "left",
                      conditionalPanel(condition = "output.dimredoutput",
                                       h4("Add new annotations to your data:")
                      ),
                      conditionalPanel(condition = "output.dimredoutput",
                                       actionButton("add_annotation", "Add new annotation!",
                                                    icon = icon("pen"))
                      )),
               column(3,align = "left",
                      conditionalPanel(condition = "output.dimredoutput",
                                       h4("Delete an annotation from your data:")
                      ),
                      conditionalPanel(condition = "output.dimredoutput",
                                       actionButton("delete_annotation", "Delete current annotation!",
                                                    icon = icon("trash"))
                                       )),
               column(3,align = "left",
                      conditionalPanel(condition = "output.dimredoutput",
                                       h4("Store your annotations for next time:")
                      ),
                      conditionalPanel(condition = "output.dimredoutput",
                                       actionButton("store_annotations_perm", "Save annotations permanently!", icon = icon("save"))
                                       )
                      )
             ),
             
             br(),
             br(),
             br(),
             hr(),
             
             
             fluidRow(
               column(4, align = "left",

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
               
               column(6,align = "left",
                      conditionalPanel(condition = "output.rename_selected == 'assigned_clusters' && output.dimredoutput",
                                       p("This panel let's you rename clusters based on previously assigned groups of cell 
                                         clusters!")),
                      
                      conditionalPanel(condition = "output.rename_selected == 'gene_expression' && output.dimredoutput",
                                       p("This panel let's you rename cell clusters based on gene expression thresholds.
                                         This can be a powerful way of classifying cells into clusters of shared expression
                                         profiles based on known markers. Currently supports up to 2 markers!")),
                      conditionalPanel(condition = "output.rename_selected == 'cell_selection' && output.dimredoutput",
                                       p("This panel lets you rename cell clusters by directly selecting a group of cells 
                                  via plotly and then assigning a cluster name to them.")
                                       )
                      )
               ),

             
             br(),
             br(),
             br(),
             hr(),
             
             
             ## Reactive panels that change based on how the user wants to rename clusters
             fluidRow(
               column(4, align = "left",
                      conditionalPanel("output.rename_selected == 'assigned_clusters'",
                                       plotOutput("dimred_plot_rename_assigned_clusters")
                                       ),
                      
                      conditionalPanel("output.rename_selected == 'gene_expression'",
                                       plotOutput("dimred_plot_rename_expression")
                                       ),
                      
                      conditionalPanel("output.rename_selected == 'cell_selection'",
                                       plotlyOutput("cell_selection_plot")
                                       ),
                      
                      conditionalPanel("output.rename_selected == 'gene_expression'",
                                       plotOutput("rename_expression_dens", height = "200px")
                                       )
                      ),
               
               column(4, align = "left",
                      uiOutput("rename_list"),
                      
                      uiOutput("gene_exp_threshold"),
                    
                      uiOutput("plot_gene_rename_button"),
                      br(),
                      br(),
                      h4(textOutput("cells_selected"))
                      ),
               
               column(4, align = "left",
                      uiOutput("cluster_anno_text"),
                      br(),
                      uiOutput("save_anno_button")
                      )
             )
    ),
    
    ## Table with marker genes
    tabPanel("Marker genes",icon = icon("highlighter"),
             
             br(),
             br(),
             
             
             fluidRow(
               column(12,
             ## Show this text only when no marker table has been uploaded!
             p("This panel lets you recalculate marker genes using the",
                      a(href="https://github.com/immunogenomics/presto", target="_blank", "Presto package"),
                "Thanks to presto, we can recalculate markers on the fly. Depending on the size of your dataset, this can take up to 10 seconds!
               Note that only the top 500 marker genes based on the AUC value per cluster annotation are shown here
               for reactivity and speed reasons!"
                )
               )
             ),
             hr(),

            br(),
            br(),

            fluidRow(
              column(4,
                     actionButton("calc_presto_markers", "Calculate marker genes!")
              ),
              column(4,
                     downloadButton("download_marker_table", "Download marker table!"))
            ),

            br(),
            br(),
            
            dataTableOutput("presto_marker_table")

            
             # conditionalPanel(condition = "output.marker_genes_table_out",
             #                  dataTableOutput("table_marker_genes")
             
             )
    
  ) # end of TabsetPanel
  ) # end of fluidPage
) # End of shinyUI
