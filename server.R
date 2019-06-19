## Load libraries
library(shiny)
library(tidyverse)
library(ggrepel)
library(viridis)
library(cowplot)
library(feather)
library(data.table)
library(plotly)
library(RColorBrewer)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
  
  ## Feather clustering file
  shinyFileChoose(input, "feather_file", 
                  roots = volumes,
                  defaultRoot = 'Home', 
                  defaultPath = 'Postdoc/Genap/data',
                  session = session)
  
  output$filepaths <- renderPrint({
    parseFilePaths(volumes, input$feather_file)
  })
  
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
    return(dimred)
  })
  
  output$dimredoutput <- reactive({
    return(dimred())
  })
  
  outputOptions(output, 'dimredoutput', suspendWhenHidden = FALSE)
  
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
  
  user_gene <- eventReactive(input$plot_gene_button ,{
    ## Check that the gene exists in the data
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
  
  ## PLot dimensional reduction plot using plotly
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
  
  ## only show tabPanel  for dimensional reduction once the plot has been created
  hideTab(inputId = "main_page", target = "Dimensional reduction")
  observeEvent(input$feather_file, {
    showTab(inputId = "main_page", target = "Dimensional reduction",
            select = TRUE)
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
                              color = "grey"),
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
  
  output$original_cell_labels <- renderUI({
    cell_classes <- unique(dimred$cell_classification)
    selectInput(inputId = "original_label",
                choices = cell_classes,
                label="Select cell cluster to relabel!" )
  })
  
  output$new_cell_labels <- renderUI({
    textInput(inputId = "new_label",
              label = "Enter new cell label:",
              value =input$original_label)
  })
  
  
  user_cell_label <- eventReactive(input$save_new_label ,{
    ## Check that the gene exists in the data
    return(input$new_label)
  })
  
  output$test_rename <-renderText({
    user_cell_label()
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
