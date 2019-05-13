## Load libraries
library(shiny)
library(tidyverse)
library(ggrepel)
library(viridis)
library(cowplot)
library(feather)
library(data.table)

## Read in data
#clust_data <- readRDS("./data/Clustering.shiny_ready.Genap2.Rds")


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  current_data <- reactive({
    
    req(input$user_gene_clustering)
    
    ## Read in data at beginning of server
    dimred <- feather::read_feather("../../test_data/clustering_df.feather",
                                    columns = c("tSNE_1","tSNE_2","classification",input$user_gene_clustering))
    
    return(dimred)
    
  })
  
  ## Calculate centroids of cluster
  centroids <- reactive({
    
    req(current_data())
    
    ## Calculate centers of each cluster
    centers <- current_data() %>%
      dplyr::group_by(classification) %>%
      summarise_at(vars(tSNE_1,tSNE_2),funs(mean(., na.rm=TRUE)))
    
    centers <- as.data.frame(centers)
    
    return(centers)
  })
   

  ## UMAP clustering with cell types
  output$tsne_plot_cluster <- renderPlot({
    
    req(current_data())
    req(centroids())
    
    tsne_plot <- ggplot(current_data(),aes(tSNE_1,tSNE_2)) +
      geom_point(aes(fill = classification),
                 size = input$point_size,
                 shape = 21,
                 alpha = 0.75,
                 colour = "black") 
    
    if(input$show_labels == TRUE){
      tsne_plot <- tsne_plot +
        geom_label_repel(data=centroids(), x = centroids()[, 2], y = centroids()[, 3], label = centroids()[,1],
                         fontface = 'bold', color = 'white', point.padding= FALSE,
                         fill = "darkgrey", size = 6, alpha = 0.8) +
        theme(legend.position = "none")
      
    }else{
      tsne_plot <- tsne_plot +
        theme(legend.position = "bottom")
    }
    
    tsne_plot
    
  })
  
  ## Plot UMAP plot with gene expression
  output$tsne_plot_gene_expression <- renderPlot({
    
    req(input$user_gene_clustering)
    
    ggplot(current_data(),aes(tSNE_1,tSNE_2)) +
      geom_point(data = subset(current_data(), get(input$user_gene_clustering) == 0 | get(input$user_gene_clustering) < 0),
                 size = input$point_size,
                 colour = "grey",
                 alpha = 0.5) +
      geom_point(data = subset(current_data(), get(input$user_gene_clustering) > 0),
                 aes(colour = get(input$user_gene_clustering)),
                 size = input$point_size,
                 alpha = 0.5) +
      scale_colour_viridis(option=input$color_palette,
                           name="Norm. Expr.") +
      labs(colour = "Norm. Expr.",
           x = "UMAP1",
           y = "UMAP1",
           title = input$gene) +
      theme(legend.title = element_text(size = 10,face="bold"),
            legend.key.size= unit(0.8, "cm"),
            legend.text = element_text(size = 10),
            axis.title.x = element_text(size =14,colour="black", face="bold"),
            axis.title.y = element_text(size =14,colour="black", face="bold"),
            strip.background=element_blank(),
            strip.text = element_blank(),
            legend.position = "right")
  })
  
  marker_table <- reactive({
    
    marker <- fread("../../test_data/marker_genes.txt")
    marker <- marker %>%
      select(-V1)
    
    marker
  })
  
  
  ## Table with marker genes to select for GeneonTSNEplot
  output$table_marker_genes <- renderDataTable(
    marker <- marker_table(),
    server=TRUE,
    caption = 'Table 1: Marker genes for all cell types',
    filter = 'top',
    selection = 'single'
  )
  
  
})
