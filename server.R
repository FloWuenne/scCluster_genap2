## Load libraries
library(shiny)
library(tidyverse)
library(ggrepel)
library(viridis)
library(feather)
library(data.table)
library(plotly)
library(RColorBrewer)
library(cowplot)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  ## Volumes for testing
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
  
  ## Volumes for GenAP2 production environment
  # volumes <- c(getVolumes()())
  
  
  #### Clustering file
  ## Feather clustering file
  shinyFileChoose(input, "feather_file", 
                  roots = volumes,
                  defaultRoot = 'Home', 
                  defaultPath = 'Postdoc/Genap/data',
                  session = session)
  
  # output$filepaths <- renderPrint({
  #   parseFilePaths(volumes, input$feather_file)
  # })
  
  ## path to the uploaded feather file
  dimred_path <- reactive({
    req(input$feather_file)
    dimred_path <- parseFilePaths(volumes, input$feather_file)$datapath
    return(dimred_path)
  })
  
  dimred <- reactive({
    req(dimred_path())
    ## Read in clustering data
    dimred <- feather::read_feather(dimred_path(),
                                    columns = c("tSNE_1","tSNE_2","cell_classification","nGene","nUMI"))
    dimred <- dimred %>%
      mutate("cell_id" = paste("cl_",rownames(dimred)))
    return(dimred)
  })
  
  ## Delete?
  # output$dimredoutput <- reactive({
  #   return(dimred())
  # })
  # 
  # outputOptions(output, 'dimredoutput', suspendWhenHidden = FALSE)
  
  ####
  ## gene names file
  shinyFileChoose(input, "gene_names", 
                  roots = volumes,
                  defaultRoot = 'Home', 
                  defaultPath = 'Postdoc/Genap/data',
                  session = session)
  
  ## path to the uploaded feather file
  gene_names_path <- reactive({
    req(input$gene_names)
    gene_names_path <- parseFilePaths(volumes, input$gene_names)$datapath
    return(gene_names_path)
  })
  
  gene_names_df <-  reactive({
    req(gene_names_path())
    gene_names_df <- fread(gene_names_path())
    return(gene_names_df)
  })
  
  #### 
  ## Marker list
  ## gene names file
  shinyFileChoose(input, "marker_genes", 
                  roots = volumes,
                  defaultRoot = 'Home', 
                  defaultPath = 'Postdoc/Genap/data',
                  session = session)
  
  ## path to the uploaded feather file
  marker_genes_path <- reactive({
    req(input$marker_genes)
    marker_genes_path <- parseFilePaths(volumes, input$marker_genes)$datapath
    return(marker_genes_path)
  })
  
  marker_genes_table <-  reactive({
    req(marker_genes_path())
    marker_genes_table <- fread(marker_genes_path())
    return(marker_genes_table)
  })
  
  ####
  ## User clustering solutions
  shinyFileChoose(input, "user_cluster_file",
                  roots = volumes,
                  defaultRoot = 'Home',
                  defaultPath = 'Postdoc/Genap/data',
                  session = session)

  ## path to the uploaded feather file
  user_cluster_path <- reactive({
    req(input$user_cluster_file)
    user_cluster_path <- parseFilePaths(volumes, input$user_cluster_file)$datapath
    return(user_cluster_path)
  })

  user_clusterings_from_file <- reactive({
    req(user_cluster_path())
    ## Read in clustering data
    user_cluster_labels <- feather::read_feather(user_cluster_path())
    return(user_cluster_labels)
  })
  
  output$cluster_columns <- renderPrint({
    req(user_clusterings_from_file())
    return(colnames(user_clusterings_from_file())[1:3])
  })

  
  ## Check that the gene exists in the data
  user_gene <- eventReactive(input$plot_gene_button ,{
    return(input$user_gene_clustering)
  })
  
  ## Loads the expression of the requested gene and adds it to the cell embeddings
  dimred_exp <- reactive({
    
    req(dimred())
    
    validate(
      need(user_gene() %in% gene_names_df()$genes,
           message = "Please enter a valid gene name")
    )
    
    gene_exp <- feather::read_feather(dimred_path(),
                                      columns = c(user_gene()))
    
    dimred_exp_df  <- cbind(dimred(),gene_exp)
    return(dimred_exp_df)
  }) 

  ## Determine color palette based number of clusters
  discrete_color_palette <- reactive({
    
    req(dimred())
    
    n_clusters <- length(unique(dimred()$cell_classification))
    
    if(n_clusters <= 9){
    colors <- brewer.pal(name = "Set1",n = n_clusters)
    }else{
    colors <- c(brewer.pal(name = "Set1",n = 9),brewer.pal(name = "Set2",n = 8))
    }
    
    return(colors)
  })
  
  ## Plot dimensional reduction plot using plotly
  output$dimred_plot_plotly <- renderPlotly({

    req(dimred())
    
    dimred_plot <- plot_ly(dimred(),
                         x = ~tSNE_1,
                         y = ~tSNE_2,
                         alpha  = 0.75,
                         type = "scattergl",
                         mode = "markers",
                         hoverinfo = 'none',
                         # text = ~paste('nGene: ', nGene, '\n',
                         #               'nUMI: ', nUMI),
                         marker = list(size = input$point_size),
                         color = ~cell_classification,
                         colors = discrete_color_palette()) 
    
    return(dimred_plot)
  })
  
  
  ## Variation of the gene expression plot using plotly
  output$dimred_gene_plot_plotly <- renderPlotly({
    
    req(dimred_exp())

    dimred_gene_plot <- plot_ly(subset(dimred_exp(),get(user_gene()) > 0),
                         x = ~tSNE_1,
                         y = ~tSNE_2,
                         alpha  = 1,
                         type = "scattergl",
                         mode = "markers",
                         marker = list(size = input$point_size),
                         color = ~get(user_gene()),
                         name = 'Gene expressed',
                         hoverinfo = 'none',
                         colors = viridis_pal(option = input$dimred_color_palette)(6)) %>%
      add_trace(p,
                data = subset(dimred_exp(),get(user_gene()) == 0),
                x = ~tSNE_1,
                y = ~tSNE_2,
                type = 'scatter', 
                mode = 'markers', 
                marker = list(size = input$point_size,
                              color = 'grey'),
                hoverinfo = 'none',
                name = 'Not expressed') %>%
      layout(title=user_gene()) %>%
      colorbar(title = "Normalized expression")
    
    
    return(dimred_gene_plot)
    
  })
  
  ## Violin plot based on plotly
  output$vlnplot_user_gene_plotly <- renderPlotly({
    
    req(dimred_exp())
    req(user_gene())
    
    vln_plot <- plot_ly(dimred_exp(),
            x = ~cell_classification,
            y = ~get(user_gene()),
            color = ~cell_classification,
            type = 'violin',
            box = list(
              visible = T),
            meanline = list(
              visible = T),
            colors = discrete_color_palette()
    )
    
    return(vln_plot)
  })
  
  ## Table with marker genes to select for GeneonTSNEplot
  output$table_marker_genes <- renderDataTable(
    {
    datatable(marker <- marker_genes_table(),
              caption = 'Table 1: Marker genes for all cell classifications',
              filter = 'top',
              selection = 'single') %>%
      formatRound(digits = c(2), columns = c(2)) %>% 
      formatStyle(columns = c(3:9), 'text-align' = 'centers')
  })
  
  ## Renaming 
  output$rename_list <- renderUI({
    
    ## If user wants to relabel assigned clusters
    if(input$rename_method == "assigned_clusters"){
      req(dimred())
      cell_classes <- unique(dimred()$cell_classification)
      selectInput(inputId = "rename_cluster_highlight",
                  choices = cell_classes,
                  label="Select cell cluster to relabel!",
                  selected = cell_classes[1])
    
      ## if user wants to label cells based on gene expression
    }else if(input$rename_method == "gene_expression"){
      req(gene_names_df())
      
      textInput("genes_renaming", label = "Enter Genesymbol", "GAPDH")
      
      }else if(input$rename_method == "cell_selection"){
        req(gene_names_df())
        req(dimred_exp())
        selectizeInput(inputId = "cells_selected", 
                       label = "Which genes would you like to use to select cells? (single-gene or comma separated list!", 
                       choices = gene_names_df()$genes[1:5], 
                       selected = NULL, multiple = FALSE,
                       options = NULL)
      }
    })
  
  output$rename_selected <- reactive({
    return(input$rename_method)
    })
  outputOptions(output, "rename_selected", suspendWhenHidden = FALSE)

  
  ## ggplot to highlight which cluster will be renamed
  output$dimred_plot_rename_assigned_clusters <- renderPlot({
    
    req(dimred())
    req(input$rename_cluster_highlight)
    
    dimred_plot <- ggplot(dimred(),aes(tSNE_1,tSNE_2)) +
      geom_point(data = subset(dimred(),cell_classification != input$rename_cluster_highlight),
                 colour = "darkgrey",
                 size = 3,
                 alpha  = 0.75) +
      geom_point(data = subset(dimred(),cell_classification == input$rename_cluster_highlight),
                 colour = "red",
                 size = 3,
                 alpha  = 1) +
      labs(x = "Dimension 1",
           y = "Dimension 2",
           title = paste("Original cluster:",input$rename_cluster_highlight))
    
    return(dimred_plot)
  })
  
  output$plot_gene_rename_button <- renderUI({
    req(gene_names_df())
    if(input$rename_method == "gene_expression"){
      actionButton("gene_rename_button", "Plot gene!")
    }
  })
  
  output$gene_exp_threshold <- renderUI({
    req(gene_names_df())
    
    if(input$rename_method == "gene_expression"){
      sliderInput(inputId = "gene_thresh_selected", 
                  label ="Select the gene expression threshold!", 
                  min = 0, max = 10,
                  value = 1, step = 0.1)
    }
  })
  
  user_genes_renaming_list <- eventReactive(input$gene_rename_button ,{
    ## Check that the gene exists in the data
    gene_vector <- unlist(strsplit(input$genes_renaming,split=","))
    return(gene_vector)
  })
  
  ## dimensional reduction with genes to plot for expression relabeling
  dimred_exp_rename <- reactive({
    req(dimred())
    req(user_genes_renaming_list())
    
    validate(
      need(user_genes_renaming_list() %in% gene_names_df()$genes,
           message = "Please enter a valid gene name")
    )
    
    gene_exp <- feather::read_feather(dimred_path(),
                                      columns = c(user_genes_renaming_list()))
    
    dimred_exp_df  <- cbind(dimred(),gene_exp)
    dimred_exp_df <- dimred_exp_df %>%
      gather("gene","expression",user_genes_renaming_list())
    
    return(dimred_exp_df)
  })
  
  ## ggplot object to highlight genes to rename clusters
  output$dimred_plot_rename_expression <- renderPlot({
    req(dimred_exp_rename())
    req(user_genes_renaming_list())
    req(input$gene_thresh_selected)
    
    dimred_plot <- ggplot(dimred_exp_rename(),aes(tSNE_1,tSNE_2)) +
      geom_point(data = subset(dimred_exp_rename(),expression < as.numeric(input$gene_thresh_selected)),
                 colour = "darkgrey",
                 size = 3,
                 alpha  = 0.75) +
      geom_point(data = subset(dimred_exp_rename(),expression >= as.numeric(input$gene_thresh_selected)),
                 colour = "red",
                 size = 3,
                 alpha  = 1) +
      labs(x = "Dimension 1",
           y = "Dimension 2",
           title = paste("Genes:",user_genes_renaming_list()))
    
    return(dimred_plot)
  })
  
  ## Print how many cells have been selected by the respective method
  output$cells_exp_selected <-  renderText({
    if(input$rename_method == "assigned_clusters"){
      req(dimred())
      cells <- subset(dimred(),cell_classification == input$rename_cluster_highlight)
      cells <- cells$cell_id
      print(paste("You have selected:",length(cells),"cells in the current assigned cluster.",sep=" "))
    }else if(input$rename_method == "gene_expression"){
      req(input$gene_thresh_selected)
      req(user_genes_renaming_list())
      cells <- subset(dimred_exp_rename(),expression >= as.numeric(input$gene_thresh_selected))
      cells <- cells$cell_id
    print(paste("You have selected:",length(cells),"cells based on the expression of ",length(user_genes_renaming_list()),"genes.",sep=" "))
    }
  }) 
  
  ## ggplot object to highlight genes to rename clusters
  output$rename_expression_hist <- renderPlot({
    req(dimred_exp_rename())
    req(user_genes_renaming_list())
    req(input$gene_thresh_selected)
    
    dimred_exp_hist <- ggplot(dimred_exp_rename(),aes(expression)) +
      geom_density(color = "darkgrey") +
      geom_vline(xintercept = as.numeric(input$gene_thresh_selected),
                 col = "red", linetype = 2)
    
    return(dimred_exp_hist)
  })
  
  output$brush <- renderPrint({
    req(output$dimred_plot_plotly)
    d <- event_data("plotly_selected")
    if (!is.null(d)) d
  })
  
  output$new_cell_labels <- renderUI({
    textInput(inputId = "new_label",
              label = "Enter new cell label:",
              value =input$original_label)
  })
  
  
  #### Controls for renaming clusters
  
  ## List of available clustering solutions
  output$choices_clusterings_rename <- renderUI({

    ## Original cluster label is cell_classification
    original_cluster_labels <- "cell_classification"
    
    ## Check if user loaded a clustering solutions file 
    if(length(colnames(user_clusterings_from_file())) > 0){
      other_clusterings <- colnames(user_clusterings_from_file())[1:3]
      all_clusterings <- c(original_cluster_labels,other_clusterings)
    }else{
      all_clusterings <- c(original_cluster_labels)
    }
    
    ## Check if the user has created other clusterings
    selectInput(inputId = "clustering_to_rename", 
                label = "Which clustering do you want to rename?",
                choices = all_clusterings)

  })
  
  clustering_solution <- eventReactive(input$start_clustering_solution ,{
    ## Create new feather file with 1 column that at the start is equivalent 
    ## to the current clustering
    current_clusters <- 
    
    return(input$user_gene_clustering)
  })
  
  

  
  
  
  
  
  
  ################### old code
  
  ## UMAP clustering with cell types
  output$tsne_plot_cluster <- renderPlot({
    
    req(dimred())
    
    tsne_plot <- ggplot(dimred(),aes(tSNE_1,tSNE_2,fill = cell_classification)) +
      geom_point(size = input$point_size,
                 shape = 21,
                 alpha = 0.75,
                 colour = "black") 
    
    if(length(unique(dimred()$cell_classification)) <= 8){
      tsne_plot <- tsne_plot + scale_fill_brewer(palette = "Set2") 
    }else if(length(unique(dimred()$cell_classification)) <= 12){
      tsne_plot <- tsne_plot + scale_fill_brewer(palette = "Set2") 
    }
    
    if(input$show_labels == TRUE){
      tsne_plot <- tsne_plot +
        # geom_label_repel(data= cluster_centers, x = cluster_centers[, 2], y = cluster_centers[, 3], label = cluster_centers[,1],
        #                  fontface = 'bold', color = 'white', point.padding= FALSE,
        #                  fill = "darkgrey", size = 6, alpha = 0.8) +
        theme(legend.position = "none")
      
    }else{
      tsne_plot <- tsne_plot +
        theme(legend.position = "bottom")
    }
    
    return(tsne_plot)
    
  })
    
  
    ## Plot cell embeddings plot with gene expression
    output$tsne_plot_gene_expression <- renderPlot({
      
      req(dimred_exp())
    
      mapping_expression <- ggplot(dimred_exp(),aes(tSNE_1,tSNE_2)) +
        geom_point(data = subset(dimred_exp(), get(user_gene()) == 0 | get(user_gene()) < 0),
                   size = input$point_size,
                   colour = "grey",
                   alpha = 0.5) +
        geom_point(data = subset(dimred_exp(), get(user_gene()) > 0),
                   aes(colour = get(user_gene())),
                   size = input$point_size,
                   alpha = 0.5) +
        scale_colour_viridis(option=input$dimred_color_palette,
                             name="Norm. Expr.") +
        labs(colour = "Norm. Expr.",
             x = "UMAP1",
             y = "UMAP1",
             title = user_gene()) +
        theme(legend.title = element_text(size = 10,face="bold"),
              legend.key.size= unit(0.8, "cm"),
              legend.text = element_text(size = 10),
              axis.title.x = element_text(size =14,colour="black", face="bold"),
              axis.title.y = element_text(size =14,colour="black", face="bold"),
              strip.background=element_blank(),
              strip.text = element_blank(),
              legend.position = "right") 
      
      return(mapping_expression)
    
      })
    
    
    output$vlnplot_user_gene <- renderPlot({
      
      req(dimred_exp())
      
      vln_plot <- ggplot(dimred_exp(),aes(cell_classification,get(user_gene()),fill=cell_classification)) +
        geom_violin() +

        labs(x="Cell cluster",
             y = "Normalized Expression",
             title = user_gene()) + 
        theme(legend.position = "none")
      
      if(length(unique(dimred_exp()$cell_classification)) <= 8){
        vln_plot <- vln_plot + scale_fill_brewer(palette = "Set2") 
      }else if(length(unique(dimred_exp()$cell_classification)) <= 12){
        vln_plot <- vln_plot + scale_fill_brewer(palette = "Set2") 
      }
      
      return(vln_plot)
    })
    
})
