## Data loading
## Default file naming convention
## 1) Clustering file: shiny_clustering_file.feather
## 2) 

## Load libraries
library(shiny)
library(tidyverse)
library(ggrepel)
library(viridis)
library(feather)
library(data.table)
library(plotly)
library(RColorBrewer)
library(shinyFiles)
library(shinyalert)
library(Matrix)
library(presto)
library(cowplot)
library(gtools)
library(DT)

theme_set(  theme_cowplot())

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  color_palette <- reactive({
    
    n_annotations <- length(unique(all_annotations()[,input$annotations_to_plot][[1]]))
    
    if(n_annotations <= 20){
    colors <- c("#3ABFEF", "#E877BB", "#C61E8E", "#EF2E25", "#7958A5", "#AC6EAE", "#CD9CC8", "#533F89",
                               "#D1B38B", "#676867", "#A5D171","#81A568", "#098040", "#3953A5", "#FFD03F", "#EF6507", "#6479F2",
                               "#438AC9", "#FFA63F", "#DB91C8")
    }else{
    colors <- "Set1"
    }
  })
  
  output$genap_logo <- renderImage({

    # Return a list containing the filename
    list(src = "./img/GenAP_powered_reg.png",
         contentType = 'image/png',
         width = "100%",
         height = "100%",
         alt = "This is alternate text")
  }, deleteFile = FALSE)
  
  ## Volumes for testing
  volumes <- c("FTP" = "/ftp",
               Home = fs::path_home())
  
  default_home <- 'FTP'
  
  ## File directory
  shinyDirChoose(input, "file_dir", 
                 roots = volumes, 
                 session = session, 
                 #defaultRoot = default_home,
                 #defaultPath = default_path,
                 restrictions = system.file(package = "base"))
  
  ## Get the location of the selected folder as a reactive variable
  file_dir_path <- reactive({
    req(input$file_dir)
    this_path <- parseDirPath(volumes, input$file_dir)
    #this_path <- parseDirPath(c("FTP" = "/ftp"), input$file_dir)
    
    ## Path for work machine (#work)
    #this_path <- parseDirPath(c("Home" = paste(fs::path_home(),sep="/")), input$file_dir)
    return(this_path)
  })
  
  output$current_dataset <- renderText({
    return(paste(tail(input$file_dir[[1]],n = 1),sep=""))
  })

  ## path to the uploaded feather file
  dimred_path <- reactive({
      req(file_dir_path())
      dimred_file <- paste(file_dir_path(),"/","shiny_clustering_file.feather",sep="")

    return(dimred_file)
  })
  
  ## Quickly check whether the user has used tSNE or UMAP for clustering
  dimred_method <- reactive({
    req(dimred_path())
    
    feather_file <- feather(dimred_path())
    first_line <- feather_file[1,]
    
    if("tSNE_1" %in% colnames(first_line)){
      method <- "tSNE"
    }else if("UMAP_1" %in% colnames(first_line)){
      method <- "UMAP"
    }
    return(method)
  })
  
  # Print method the user was using
  output$dimred_met <- renderPrint({
    dimred_method()
  })
  
  ## path to the uploaded feather file
  presto_path <- reactive({
    req(file_dir_path())
    presto_file <- paste(file_dir_path(),"/","shiny_clustering.sparse_presto.rds",sep="")
    return(presto_file)
  })
  
  dimred <- reactive({
    req(dimred_path())
    req(input$annotations_to_plot)
    req(dimred_method())
    
    ## Need to add support for UMAP and split binding of a  nnotations to dimred for only reading feather standard file
    ## once when annotations change.
    ## Idea for UMAP: quickly read file that contains parameters?! 
    
    ## Read in clustering data
    if(dimred_method() == "tSNE"){
      dimred <- feather::read_feather(dimred_path(),
                                      columns = c("tSNE_1","tSNE_2","cell_id"))
    }else if(dimred_method() == "UMAP"){
      dimred <- feather::read_feather(dimred_path(),
                                      columns = c("UMAP_1","UMAP_2","cell_id"))
    }


    selected_annotation <- all_annotations()[,input$annotations_to_plot]
    dimred <- cbind(dimred,selected_annotation)
    
    return(dimred)
  })
  
  ## Reactive output for enabling conditionalPanel based on file input
  output$dimredoutput <- reactive({
    return(dimred())
  })

  outputOptions(output, 'dimredoutput', suspendWhenHidden = FALSE)
  
  ####
  ## gene names file
  shinyFileChoose(input, "gene_names", 
                  roots = volumes,
                  defaultRoot = default_home, 
                  defaultPath = default_path,
                  session = session)
  
  ## path to the uploaded feather file
  gene_names_path <- reactive({

    req(file_dir_path())
    gene_names_path <- paste(file_dir_path(),"/","shiny_gene_names.tsv",sep="")

    return(gene_names_path)
  })
  
  gene_names_df <-  reactive({
    req(gene_names_path())
    gene_names_df <- fread(gene_names_path())
    return(gene_names_df)
  })

  ####
  ## User clustering solutions
  shinyFileChoose(input, "user_cluster_file",
                  roots = volumes,
                  defaultRoot = default_home,
                  defaultPath = default_path,
                  session = session)

  ## path to the uploaded feather file
  user_cluster_path <- reactive({
    req(file_dir_path())
    user_cluster_path <- paste(file_dir_path(),"/","shiny_user_clustering.feather",sep="")
    return(user_cluster_path)
  })

  user_clusterings_from_file <- reactive({
    req(user_cluster_path())
    ## Read in clustering data
    user_cluster_labels <- feather::read_feather(user_cluster_path())
    return(user_cluster_labels)
  })
  
  ## Reactive output for enabling conditionalPanel based on file input
  output$user_cluster_labels <- reactive({
    return(user_clusterings_from_file())
  })
  
  outputOptions(output, 'user_cluster_labels', suspendWhenHidden = FALSE)
  
  
  ## Check that the gene exists in the data
  user_gene <- eventReactive(input$plot_gene_button,{
    return(input$user_gene_clustering)
  })
  
  user_gene_alert <- observeEvent(input$plot_gene_button, {
    if(is.null(all_annotations())){
      shinyalert("Error!", "Please upload a dataset first!", type = "error")
    }
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

  
  ## Plot dimensional reduction plot using plotly
  output$dimred_plot_plotly <- renderPlotly({

    req(dimred())
    req(input$annotations_to_plot)
    req(color_palette())
    
    f <- list(
      size = 18,
      color = "black"
    )
    
    x <- list(
      title = paste(dimred_method()," 1",sep=""),
      titlefont = f
    )
    y <- list(
      title = paste(dimred_method()," 2",sep=""),
      titlefont = f
    )
    
    dimred_plot <- plot_ly(dimred(),
                         # x = ~tSNE_1,
                         # y = ~tSNE_2,
                         x = ~get(paste(dimred_method(),"_1",sep="")),
                         y = ~get(paste(dimred_method(),"_2",sep="")),
                         alpha  = 0.75,
                         type = "scattergl",
                         mode = "markers",
                         hoverinfo = 'text',
                         text = ~get(input$annotations_to_plot),
                         marker = list(size = input$point_size),
                         color = ~get(input$annotations_to_plot),
                         colors = color_palette()) %>%
      layout(xaxis = x, yaxis = y)
    
    return(dimred_plot)
  })
  
  
  ## Variation of the gene expression plot using plotly
  output$dimred_gene_plot_plotly <- renderPlotly({
    
    req(dimred_exp())
    
    f <- list(
      size = 18,
      color = "black"
    )
    
    x <- list(
      title = paste(dimred_method()," 1",sep=""),
      titlefont = f
    )
    y <- list(
      title = paste(dimred_method()," 2",sep=""),
      titlefont = f
    )

    dimred_gene_plot <- plot_ly(subset(dimred_exp(),get(user_gene()) > 0),
                         # x = ~tSNE_1,
                         # y = ~tSNE_2,
                         x = ~get(paste(dimred_method(),"_1",sep="")),
                         y = ~get(paste(dimred_method(),"_2",sep="")),
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
                # x = ~tSNE_1,
                # y = ~tSNE_2,
                x = ~get(paste(dimred_method(),"_1",sep="")),
                y = ~get(paste(dimred_method(),"_2",sep="")),
                type = 'scatter', 
                mode = 'markers', 
                marker = list(size = input$point_size,
                              color = 'grey'),
                hoverinfo = 'none',
                name = 'Not expressed') %>%
      layout(title=user_gene(),
             xaxis = x, yaxis = y) %>%
      colorbar(title = "Normalized expression")
    
    
    return(dimred_gene_plot)
    
  })
  
  ## Violin plot based on plotly
  output$vlnplot_user_gene_plotly <- renderPlotly({
    
    req(dimred_exp())
    req(user_gene())
    req(input$annotations_to_plot)
    
    x <- list(
      title = "Cluster"
    )
    y <- list(
      title = "Normalized expression"
    )
    
    vln_plot <- plot_ly(dimred_exp(),
            x = ~get(input$annotations_to_plot),
            y = ~get(user_gene()),
            color = ~get(input$annotations_to_plot),
            type = 'violin',
            box = list(
              visible = T),
            meanline = list(
              visible = T)) %>%
      layout(xaxis = x, yaxis = y)
    
     
    return(vln_plot)
  })
  
  presto_alert <- observeEvent(input$calc_presto_markers, {
    if(is.null(all_annotations())){
      shinyalert("Error!", "Please upload a dataset first!", type = "error")
    }
  })
  
  presto_marker_genes <- eventReactive(input$calc_presto_markers,{
    
    req(presto_path())
    req(input$annotations_to_plot)
    req(all_annotations())
    
    withProgress(message = 'Calculating Markers...', value = 0, {
      
      ## Load complete dataset for marker calculation
      #dimred <- feather::read_feather(dimred_path())
      dimred_genes <- readRDS(presto_path())
      
      # Increment the progress bar, and update the detail text.
      incProgress(0.5, detail = paste("Done reading full expression matrix..."))
      
      ## Subset complete matrix to only contain genes
      #dimred_genes <- dimred[,gene_names_df()$genes]
      annotations <- all_annotations()
      selected_annotation <- annotations[,input$annotations_to_plot]
      
      ## Run presto (if annotation length matches matrix size)
      if(ncol(dimred_genes) == length(selected_annotation[[1]])){
        #presto_results <- wilcoxauc(as(t(dimred_genes), "sparseMatrix"),selected_annotation[[1]])
        presto_results <- wilcoxauc(dimred_genes,selected_annotation[[1]])
        presto_results <- presto_results %>%
          mutate("pct_diff" = pct_in - pct_out) %>%
          subset(auc >= 0.5 & avgExpr > 0 & pct_diff > 0) %>%
          top_n(500,wt = "auc") 
        
        presto_results$group <- as.factor(presto_results$group)

        presto_results <- presto_results %>%
          group_by(group) %>%
          arrange(desc(logFC))
        
        # Increment the progress bar, and update the detail text.
        incProgress(0.8, detail = paste("Presto run finished!"))
        
      }else{
        row <- paste("Warning, matrix has",nrow(dimred_genes),"and annotations have:",
                     length(selected_annotation),"mismatch!",sep="")
        
        presto_results <- data.frame("warning" = row,
                                     "anno_head" = selected_annotation[[1]])
      }
    
    
    # Increment the progress bar, and update the detail text.
    incProgress(1, detail = paste("Marker calculation done!"))
    
    })

    return(presto_results)
    
  })
  
  output$presto_marker_table <- renderDataTable({
    DT::datatable(marker <- presto_marker_genes(),
              caption = 'Table 1: Presto Marker genes for selected annotation',
              filter = 'top',
              selection = 'single',
              rownames= FALSE) %>%
      formatRound(digits = c(2), columns = c(3:11)) %>%
      formatStyle(columns = c(1:11), 'text-align' = 'centers')
  })
  
  ## old table for marker genes precalculated from Seurat or Scanpy
  # ## Table with marker genes to select for GeneonTSNEplot
  # output$table_marker_genes <- renderDataTable({
  #     
  #   ## expression_matrix
  #     
  #   # Precomputed marker genes table
  #   datatable(marker <- marker_genes_table(),
  #             caption = 'Table 1: Marker genes for all cell classifications',
  #             filter = 'top',
  #             selection = 'single') %>%
  #     formatRound(digits = c(2), columns = c(2)) %>%
  #     formatStyle(columns = c(3:9), 'text-align' = 'centers')
  # })
  
  ## Gene to plot in Modify & add annotations view 
  gene_names_selection <- eventReactive(input$rename_method == "gene_expression",{
    req(gene_names_df())
    
    gene_names_selection <- unique(gene_names_df()$genes)
    
    updateSelectizeInput(session, 'gene_renaming', 
                         choices = c("Type your gene"='',gene_names_selection), 
                         server = TRUE,
                         selected = NULL)
  })
  
  ## Gene to plot in Dimensional reduction view
  output$user_gene_plot_list <- renderUI({
    req(dimred_gene_sel())
    selectizeInput(
      inputId = 'user_gene_clustering', 
      label = 'Select the genes to plot', 
      choices = NULL , multiple = FALSE)
  })
  
  dimred_gene_sel <- reactive({

    req(gene_names_df())
    gene_names_available <- unique(gene_names_df()$genes)
    
    updateSelectizeInput(session, 'user_gene_clustering', 
                         choices = c("Type your gene"='',gene_names_available), 
                         server = TRUE,
                         selected = NULL)
  })
  


  ## Dimensional reduction
  
  
  ## Renaming 
  output$rename_list <- renderUI({
    
    ## If user wants to relabel assigned clusters
    if(input$rename_method == "assigned_clusters"){
      req(all_annotations())
      annotations <- all_annotations()
      cell_classes <- mixedsort(unique(annotations[,input$annotations_to_plot][[1]]))
      selectInput(inputId = "rename_cluster_highlight",
                  choices = cell_classes,
                  label="Select cell cluster to relabel!")
    
      ## if user wants to label cells based on gene expression
    }else if(input$rename_method == "gene_expression"){
      req(gene_names_selection())
      
      selectizeInput(
        inputId = 'gene_renaming', 
        label = 'Select the genes to plot', 
        choices = NULL , multiple = FALSE)
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
    
    dimred_plot <- ggplot(dimred(),
                          aes(get(paste(dimred_method(),"_1",sep="")),
                              get(paste(dimred_method(),"_2",sep="")))) +
      geom_point(data = subset(dimred(),get(input$annotations_to_plot) != input$rename_cluster_highlight),
                 colour = "darkgrey",
                 size = 3,
                 alpha  = 0.75) +
      geom_point(data = subset(dimred(),get(input$annotations_to_plot) == input$rename_cluster_highlight),
                 colour = "red",
                 size = 3,
                 alpha  = 1) +
      labs(x = paste(dimred_method()," 1",sep=""),
           y = paste(dimred_method()," 2",sep=""),
           title = paste("Cluster:",input$rename_cluster_highlight))
    
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
    req(dimred_exp_rename())
    dimred_exp_rename <- subset(dimred_exp_rename(),expression > 0)
    
    if(input$rename_method == "gene_expression"){
      sliderInput(inputId = "gene_thresh_selected", 
                  label ="Select the gene expression threshold!", 
                  min = round(min(dimred_exp_rename$expression),2), max = round(max(dimred_exp_rename$expression),2),
                  value = 1, step = 0.1)
    }
  })
  
  
  ## dimensional reduction with genes to plot for expression relabeling
  dimred_exp_rename <- eventReactive(req(input$gene_rename_button),{
    req(dimred())
    req(input$gene_renaming)
    req(input$gene_rename_button)
    
    validate(
      need(input$gene_renaming %in% gene_names_df()$genes,
           message = "Please enter a valid gene name")
    )
    
    gene_exp <- feather::read_feather(dimred_path(),
                                      columns = c(input$gene_renaming))
    
    dimred_exp_df  <- cbind(dimred(),gene_exp)
    dimred_exp_df <- dimred_exp_df %>%
      gather("gene","expression",input$gene_renaming)
    
    return(dimred_exp_df)
  })
  
  ## ggplot object to highlight genes to rename clusters
  output$dimred_plot_rename_expression <- renderPlot({
    req(dimred_exp_rename())
    req(input$gene_renaming)
    req(input$gene_thresh_selected)
  
    
    dimred_plot <- ggplot(dimred_exp_rename(),aes(get(paste(dimred_method(),"_1",sep="")),
                                                  get(paste(dimred_method(),"_2",sep="")))) +
      geom_point(data = subset(dimred_exp_rename(),expression < as.numeric(input$gene_thresh_selected)),
                 colour = "darkgrey",
                 size = 3,
                 alpha  = 0.5) +
      geom_point(data = subset(dimred_exp_rename(),expression >= as.numeric(input$gene_thresh_selected)),
                 colour = "red",
                 size = 3,
                 alpha  = 0.75) +
      labs(x = paste(dimred_method()," 1",sep=""),
           y = paste(dimred_method()," 2",sep=""),
           title = paste("Gene:",input$gene_renaming),sep=" ")
    
    return(dimred_plot)
  })

  ## Print how many cells have been selected by the respective method
  output$cells_selected <-  renderText({
    if(input$rename_method == "assigned_clusters"){
      
      req(dimred())
      req(input$rename_cluster_highlight)
      
      cells <- subset(dimred(),get(input$annotations_to_plot) == input$rename_cluster_highlight)
      cells <- cells$cell_id
      print(paste("You have selected:",length(cells),"cells in the current assigned cluster.",sep=" "))
      
    }else if(input$rename_method == "gene_expression"){
      
      req(input$gene_renaming)
      
      all_cells_expressing  <- subset(dimred_exp_rename(),
                                      expression >= input$gene_thresh_selected)
      
      all_cells_expressing <- unique(all_cells_expressing$cell_id)
      print(paste("You have selected:",length(all_cells_expressing),"cells based on the expression of ",input$gene_renaming,sep=" "))
    } else if(input$rename_method == "cell_selection"){
      req(event_data("plotly_selected"))
      cells_selected <- event_data("plotly_selected")$key
      
      print(paste("You have selected:",length(cells_selected),"cells!",sep=" "))
    }
  }) 
  
  ## ggplot object to highlight genes to rename clusters
  output$rename_expression_dens <- renderPlot({
    req(dimred_exp_rename())
    req(input$gene_renaming)
    req(input$gene_thresh_selected)
    
    dimred_exp <- dimred_exp_rename()
    dimred_exp <- subset(dimred_exp,expression > 0)
    
    dimred_exp_hist <- ggplot(dimred_exp,aes(expression)) +
      geom_density(fill = "red", alpha = 1) +
      geom_vline(xintercept = as.numeric(input$gene_thresh_selected),
                 col = "black", linetype = 2,size = 1.5) +
      scale_x_continuous(breaks = round(seq(min(dimred_exp$expression), max(dimred_exp$expression), by = 1),1)) +
      labs(x = "Normalized expression",
           y = "Density",
           title= "Expression over all cells",
           subtitle = "Cells with no expression excluded") 
    
    return(dimred_exp_hist)
  })
  
  output$new_cell_labels <- renderUI({
    textInput(inputId = "new_label",
              label = "Enter new cell label:",
              value =input$original_label)
  })
  
  
  #### Controls for renaming clusters
  
  ## Create reactive value to store annotations
  all_annotations <- reactiveVal()
  
  ## Variable to hold the last selected annotation
  last_selected_annotation <- reactiveVal()
  
  ## Add annotations once file with initial annotations gets uploaded
  observeEvent(user_clusterings_from_file(),{
    all_annotations(user_clusterings_from_file())
    last_selected_annotation(colnames(user_clusterings_from_file())[1])
  })
  
  ## Add a new annotation whenever the user wants. The clusterings in this annotation will be based on the
  ## currently selected annotation
  observeEvent(input$add_annotation,{
    last_selected_annotation(input$annotations_to_plot)
    
    current_annotations <- all_annotations()
    
    ## Make new column name based on user input
    clustering_name <- input$user_added_cluster
    if(clustering_name == ""){
      ## If user doesn't enter name, generate new name automatically
      clustering_name <- paste("user_annotation",ncol(current_annotations)+1,sep="4")
    }else{
      clustering_name <- clustering_name
    }
    
    ## Add new column to annotations
    current_annotations[,clustering_name] <- current_annotations[,input$annotations_to_plot][[1]]
    all_annotations(current_annotations)
  })
  
  ## Delete currently selected annotation
  observeEvent(input$delete_annotation,{
    current_annotations <- all_annotations()
    column_names <- colnames(current_annotations)
    if(input$annotations_to_plot != "cell_classification"){
      remaining_column <- setdiff(column_names,input$annotations_to_plot)
      current_annotations <- current_annotations[,remaining_column]
    }else{
      shinyalert("STOP!", "You cannot delete the initial annotations!", type = "error")
    }
    all_annotations(current_annotations)
  })
  
  ## Rename clusters in a user selected annotation. The way this is done depends on which methods the user
  ## chooses
  observeEvent(input$save_new_annotation,{
    last_selected_annotation(input$annotations_to_plot)
    
    ## If user is renaming clusters, substitute the existing cluster name with a new one
    if(input$rename_method == "assigned_clusters"){
      current_annotations <- all_annotations()
      
      cells_to_rename <- subset(current_annotations,get(input$annotations_to_plot) == input$rename_cluster_highlight)
      
      current_annotations <- current_annotations %>%
        mutate("new_anno" = get(input$annotations_to_plot)) %>%
        mutate("new_anno" = if_else(get(input$annotations_to_plot) == input$rename_cluster_highlight,
                                    input$new_cluster_annotation,
                                    as.character(new_anno)))
      
      ## Keep order of clusters numerical and then alphabetical!
      current_annotations$new_anno <- factor(current_annotations$new_anno, mixedsort(unique(current_annotations$new_anno)))
      
      ## Reput clusterings into the original naming slot
      current_annotations[,input$annotations_to_plot] <- current_annotations[,"new_anno"]
      current_annotations$new_anno <- NULL
      
      # current_annotations[,] <- gsub(input$rename_cluster_highlight,
      #                                              input$new_cluster_annotation,
      #                                              current_annotations[,current_column][[1]])
      all_annotations(current_annotations)
      
    ## If user is renaming cells based on gene expression, mutate using if else
    } else if(input$rename_method == "gene_expression"){
      req(dimred_exp_rename())
      ## Get cells that pass the gene expression threshold
      cells_to_rename <- dimred_exp_rename() %>%
        subset(expression >= as.numeric(input$gene_thresh_selected))
      
      cells_to_rename <- cells_to_rename$cell_id
      
      current_annotations <- all_annotations()

      current_annotations <- current_annotations %>%
        mutate("new_anno" = get(input$annotations_to_plot)) %>%
        mutate("new_anno" = if_else(cell_id %in% cells_to_rename,
                                input$new_cluster_annotation,
                                as.character(new_anno)))
      
      current_annotations$new_anno <- factor(current_annotations$new_anno, mixedsort(unique(current_annotations$new_anno)))
      current_annotations[,input$annotations_to_plot] <- current_annotations[,"new_anno"]
      current_annotations$new_anno <- NULL
      all_annotations(current_annotations)
      
    } else if(input$rename_method == "cell_selection"){
      req(dimred())
      current_annotations <- all_annotations()
      cells_selected <- event_data("plotly_selected")$key
      
      current_annotations <- current_annotations %>%
        mutate("new_anno" = get(input$annotations_to_plot)) %>%
        mutate("new_anno" = if_else(cell_id %in% cells_selected,
                                    input$new_cluster_annotation,
                                    as.character(new_anno)))
      
      current_annotations$new_anno <- factor(current_annotations$new_anno, mixedsort(unique(current_annotations$new_anno)))
      current_annotations[,input$annotations_to_plot] <- current_annotations[,"new_anno"]
      current_annotations$new_anno <- NULL
      all_annotations(current_annotations)
    }

  })

  ## List of available clustering solutions
  output$available_cluster_labels <- renderUI({
    
    req(all_annotations())
    req(last_selected_annotation)

      ## Check if the user has created other clusterings
      selectInput(inputId = "annotations_to_plot", 
                  label = "Which annotations do you want to use?",
                  choices = colnames(all_annotations()[,-1]),
                  selected = last_selected_annotation()[[1]])
    
  })
  
  output$selected_annotation <- renderText({
    req(input$annotations_to_plot)
    
    return(input$annotations_to_plot)
    
  })
  
  ## List of available clustering solutions
  output$choices_clusterings_rename <- renderUI({

    req(all_annotations())

    ## Check if the user has created other clusterings
    selectInput(inputId = "clustering_to_rename", 
                label = "Which clustering do you want to rename?",
                choices = colnames(all_annotations()))

  })
  
  output$cluster_anno_text <- renderUI({
    req(dimred())
    textInput("new_cluster_annotation", label = "Enter new cluster annotation!", 10)
  })
  
  
  output$save_anno_button <- renderUI({
    req(dimred())
    actionButton("save_new_annotation", "Save new annotation label!")
  })
  
  observeEvent(input$store_annotations_perm,{
    write_feather(all_annotations(),
                    path <- paste(file_dir_path(),"/","shiny_user_clustering.feather",sep=""))
    # Show a modal when the button is pressed
    shinyalert("Success!", "Your annotations were saved", type = "success")
  })

  output$cell_selection_plot <- renderPlotly({
    
    req(dimred())
    req(input$annotations_to_plot)
    
    f <- list(
      size = 18,
      color = "black"
    )
    
    x <- list(
      title = paste(dimred_method()," 1",sep=""),
      titlefont = f
    )
    y <- list(
      title = paste(dimred_method()," 2",sep=""),
      titlefont = f
    )
  
    cell_ids <- dimred()$cell_id
    dimred_plot <- plot_ly(dimred(),
                           # x = ~tSNE_1,
                           # y = ~tSNE_2,
                           x = ~get(paste(dimred_method(),"_1",sep="")),
                           y = ~get(paste(dimred_method(),"_2",sep="")),
                           key = ~cell_ids,
                           alpha  = 0.75,
                           type = "scattergl",
                           mode = "markers",
                           hoverinfo = 'text',
                           text = ~get(input$annotations_to_plot),
                           marker = list(size = input$point_size),
                           color = ~get(input$annotations_to_plot),
                           colors = color_palette()) %>%
      layout(dragmode = "select",
             xaxis = x, yaxis = y)
    
    dimred_plot
    
  })
  
  output$brush <- renderPrint({
    d <- event_data("plotly_selected")
    if (is.null(d)) "Click and drag events (i.e., select/lasso) appear here (double-click to clear)" else d
  })
  
  # Downloadable csv of selected dataset ----
  output$download_marker_table <- downloadHandler(
    filename = function() {
      paste(input$annotations_to_plot,".tsv", sep = "")
    },
    content = function(file) {
      write.table(presto_marker_genes(), 
                  file,
                  row.names = FALSE,
                  col.names = TRUE,
                  sep="\t",
                  quote = FALSE)
    }
  )
  
  ## List 1 of annotations
  output$comp_anno_list1 <- renderUI({
    req(all_annotations())
    
    selectInput(
      inputId = 'comp_anno_1', 
      label = 'Select the genes to plot', 
      choices = colnames(all_annotations()[,-1]), multiple = FALSE)
  })
  
  ## List 2 of annotations
  output$comp_anno_list2 <- renderUI({
    req(all_annotations())
    req(input$comp_anno_1)
    
    selectInput(
      inputId = 'comp_anno_2', 
      label = 'Select the genes to plot', 
      choices = colnames(all_annotations()[,-1]), multiple = FALSE)
  })
  
  sankey_comp <- eventReactive(input$compare_annos, {
    req(all_annotations())
    req(input$comp_anno_1)
    req(input$comp_anno_2)

    annos_to_compare <- all_annotations()[,c("cell_id",input$comp_anno_1,input$comp_anno_2)]

    annos_to_compare <- all_annotations()[,c(input$comp_anno_1,input$comp_anno_2)]
    
    annos_to_compare_stats <- annos_to_compare %>% 
      group_by(get(input$comp_anno_1),get(input$comp_anno_2)) %>%
      tally() %>%
      ungroup()
    
    colnames(annos_to_compare_stats) <- c("anno1","anno2","n")
    
    annos_to_compare_stats <- annos_to_compare_stats %>%
      mutate("anno1" = paste(anno1,input$comp_anno_1,sep="_")) %>%
      mutate("anno2" = paste(anno2,input$comp_anno_2,sep="_"))
    
    joined_annos <- c(annos_to_compare_stats$anno1,annos_to_compare_stats$anno2)
    joined_annos <- unique(joined_annos)
    
    annos_to_compare_stats$IDsource=match(annos_to_compare_stats$anno1, joined_annos)-1 
    annos_to_compare_stats$IDtarget=match(annos_to_compare_stats$anno2, joined_annos)-1

    return(annos_to_compare_stats)
    })
  


  ## Sankey diagram to compare two annotations
  output$sankey_diagram <- renderPlotly({
    req(sankey_comp())
    
    joined_annos <- c(sankey_comp()$anno1,sankey_comp()$anno2)
    joined_annos <- unique(joined_annos)
    
    p <- plot_ly(
      type = "sankey",
      orientation = "h",
      
      node = list(
        label = joined_annos,
        pad = 15,
        thickness = 20,
        line = list(
          color = "black",
          width = 0.5
        )
      ),
      
      link = list(
        source = sankey_comp()$IDsource,
        target = sankey_comp()$IDtarget,
        value =  sankey_comp()$n
      )
    )



    }
  )

  
  
  
  # 
  # ################### old code
  # 
  # ## UMAP clustering with cell types
  # output$tsne_plot_cluster <- renderPlot({
  #   
  #   req(dimred())
  #   
  #   tsne_plot <- ggplot(dimred(),aes(tSNE_1,tSNE_2,fill = get(input$annotations_to_plot))) +
  #     geom_point(size = input$point_size,
  #                shape = 21,
  #                alpha = 0.75,
  #                colour = "black") 
  #   
  #   if(length(unique(dimred()[,input$annotations_to_plot])) <= 8){
  #     tsne_plot <- tsne_plot + scale_fill_brewer(palette = "Set2") 
  #   }else if(length(unique(dimred()[,input$annotations_to_plot])) <= 12){
  #     tsne_plot <- tsne_plot + scale_fill_brewer(palette = "Set2") 
  #   }
  #   
  #   if(input$show_labels == TRUE){
  #     tsne_plot <- tsne_plot +
  #       # geom_label_repel(data= cluster_centers, x = cluster_centers[, 2], y = cluster_centers[, 3], label = cluster_centers[,1],
  #       #                  fontface = 'bold', color = 'white', point.padding= FALSE,
  #       #                  fill = "darkgrey", size = 6, alpha = 0.8) +
  #       theme(legend.position = "none")
  #     
  #   }else{
  #     tsne_plot <- tsne_plot +
  #       theme(legend.position = "bottom")
  #   }
  #   
  #   return(tsne_plot)
  #   
  # })
  #   
  # 
  #   ## Plot cell embeddings plot with gene expression
  #   output$tsne_plot_gene_expression <- renderPlot({
  #     
  #     req(dimred_exp())
  #   
  #     mapping_expression <- ggplot(dimred_exp(),aes(tSNE_1,tSNE_2)) +
  #       geom_point(data = subset(dimred_exp(), get(user_gene()) == 0 | get(user_gene()) < 0),
  #                  size = input$point_size,
  #                  colour = "grey",
  #                  alpha = 0.5) +
  #       geom_point(data = subset(dimred_exp(), get(user_gene()) > 0),
  #                  aes(colour = get(user_gene())),
  #                  size = input$point_size,
  #                  alpha = 0.5) +
  #       scale_colour_viridis(option=input$dimred_color_palette,
  #                            name="Norm. Expr.") +
  #       labs(colour = "Norm. Expr.",
  #            x = "UMAP1",
  #            y = "UMAP1",
  #            title = user_gene()) +
  #       theme(legend.title = element_text(size = 10,face="bold"),
  #             legend.key.size= unit(0.8, "cm"),
  #             legend.text = element_text(size = 10),
  #             axis.title.x = element_text(size =14,colour="black", face="bold"),
  #             axis.title.y = element_text(size =14,colour="black", face="bold"),
  #             strip.background=element_blank(),
  #             strip.text = element_blank(),
  #             legend.position = "right") 
  #     
  #     return(mapping_expression)
  #   
  #     })
  #   
  #   
  #   output$vlnplot_user_gene <- renderPlot({
  #     
  #     req(dimred_exp())
  #     
  #     vln_plot <- ggplot(dimred_exp(),aes(get(input$annotations_to_plot),get(user_gene()),fill=get(input$annotations_to_plot))) +
  #       geom_violin() +
  # 
  #       labs(x="Cell cluster",
  #            y = "Normalized Expression",
  #            title = user_gene()) + 
  #       theme(legend.position = "none")
  #     
  #     if(length(unique(dimred_exp()[,input$annotations_to_plot])) <= 8){
  #       vln_plot <- vln_plot + scale_fill_brewer(palette = "Set2") 
  #     }else if(length(unique(dimred_exp()[,input$annotations_to_plot])) <= 12){
  #       vln_plot <- vln_plot + scale_fill_brewer(palette = "Set2") 
  #     }
  #     
  #     return(vln_plot)
  #   })
  #   
})
