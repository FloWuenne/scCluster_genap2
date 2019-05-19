## Load libraries
library(shiny)
library(tidyverse)
library(ggrepel)
library(viridis)
library(cowplot)
library(feather)
library(data.table)

## Read in data
## Testing locally
data_dir <- "/Users/florian_wuennemann/Postdoc/Genap/data/"
## Production Genap2 environment
#data_dir <- ""

## Read in clustering data
dimred <- feather::read_feather(paste(data_dir,"clustering_shiny.feather",sep="/"))
#columns = c("tSNE_1","tSNE_2","cell_classification",input$user_gene_clustering)

## Calculate centers of cluster labels
cluster_centers <- feather::read_feather(paste(data_dir,"cluster_info_shiny.feather",sep="/"))
cluster_centers <- as.data.frame(cluster_centers)

## Marker Table
# marker <- fread("../../test_data/marker_genes.txt")
# marker <- marker %>%
#   select(-V1)


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  observeEvent(input[["keyPressed"]], {
    input$user_gene_clustering
  })
  
  ## UMAP clustering with cell types
  output$tsne_plot_cluster <- renderPlot({
    
    tsne_plot <- ggplot(dimred,aes(tSNE_1,tSNE_2)) +
      geom_point(aes(fill = cell_classification),
                 size = input$point_size,
                 shape = 21,
                 alpha = 0.75,
                 colour = "black") 
    
    if(input$show_labels == TRUE){
      tsne_plot <- tsne_plot +
        geom_label_repel(data= cluster_centers, x = cluster_centers[, 2], y = cluster_centers[, 3], label = cluster_centers[,1],
                         fontface = 'bold', color = 'white', point.padding= FALSE,
                         fill = "darkgrey", size = 6, alpha = 0.8) +
        theme(legend.position = "none")
      
    }else{
      tsne_plot <- tsne_plot +
        theme(legend.position = "bottom")
    }
    
    return(tsne_plot)
    
  })
  
  
  observeEvent(input$plot_gene_button, {
    ## Plot cell embeddings plot with gene expression
    output$tsne_plot_gene_expression <- renderPlot({
      
      req(input$user_gene_clustering)
      
      ggplot(dimred,aes(tSNE_1,tSNE_2)) +
        geom_point(data = subset(dimred, get(input$user_gene_clustering) == 0 | get(input$user_gene_clustering) < 0),
                   size = input$point_size,
                   colour = "grey",
                   alpha = 0.5) +
        geom_point(data = subset(dimred, get(input$user_gene_clustering) > 0),
                   aes(colour = get(input$user_gene_clustering)),
                   size = input$point_size,
                   alpha = 0.5) +
        scale_colour_viridis(option=input$color_palette,
                             name="Norm. Expr.") +
        labs(colour = "Norm. Expr.",
             x = "UMAP1",
             y = "UMAP1",
             title = input$user_gene_clustering) +
        theme(legend.title = element_text(size = 10,face="bold"),
              legend.key.size= unit(0.8, "cm"),
              legend.text = element_text(size = 10),
              axis.title.x = element_text(size =14,colour="black", face="bold"),
              axis.title.y = element_text(size =14,colour="black", face="bold"),
              strip.background=element_blank(),
              strip.text = element_blank(),
              legend.position = "right")
    })
  })
  
  
  # ## Table with marker genes to select for GeneonTSNEplot
  # output$table_marker_genes <- renderDataTable(
  #   marker <- marker_table(),
  #   server=TRUE,
  #   caption = 'Table 1: Marker genes for all cell types',
  #   filter = 'top',
  #   selection = 'single'
  # )
  
  
})
