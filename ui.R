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
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      textInput("user_gene_clustering", label = "Enter Genesymbol", "Gapdh"),
      
      actionButton("plot_gene_button", "Plot gene!"),
      
      br(),
      br(),
      
      ## Dropdown menu for color panel
      selectInput(
        inputId="color_palette", label="Choose a color palette for expression!", 
        choices=c("A","B","C","D","E"),
        selected="D", multiple=FALSE, selectize=FALSE),
      
      ## Scale menu for point size 
      sliderInput(
        inputId="point_size", label="Change cell point size!", 
        min = 2, max = 10, value = 6, step = 2, round = TRUE) 
      
      # h4("Do you want to show cell type labels?"),
      # checkboxInput("show_labels",
      #               label = "Show labels",
      #               value = TRUE),
      
      # h4("Do you want to regroup or relabel cells?"),
      # checkboxInput("relabel_cells",
      #               label = "Relabel cells?",
      #               value = FALSE),
      
      # conditionalPanel(
      #   condition = "input.relabel_cells == true",
      #   uiOutput("original_cell_labels")
      # ),
      
      # conditionalPanel(
      #   condition = "input.relabel_cells == true",
      #   uiOutput("new_cell_labels"),
      #   actionButton("save_new_label", "Save new label!")
      # ),
      
      # textOutput("test_rename")

    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("tSNE Clustering",

                 plotlyOutput("tsne_plot_cluster_plotly"),
                 br(),
                 plotlyOutput("tsne_plot_gene_expression_plotly"),
                 br(),
                 plotlyOutput("vlnplot_user_gene_plotly")
                 
                 ## old plots based on ggplot
                 #withSpinner(plotOutput("tsne_plot_cluster")),
                 #plotOutput("tsne_plot_gene_expression"),
                 #plotOutput("vlnplot_user_gene")

                 
      ),
      
      tabPanel("Marker_genes",
               dataTableOutput("table_marker_genes")
      )
    
  )
  
  )
)
))
