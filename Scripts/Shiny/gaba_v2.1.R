# testing it out on the gaba dataset

library(shiny)
library(ggplot2)
library(Seurat)
library(umap)
library(dplyr)
library(patchwork)
library(stringr)
library(bslib)
library(fontawesome)
library(magick)
library(shinydashboard)
library(tidyverse)

setwd("~/Desktop/Analysis/Shiny/Test")

load("~/Desktop/Analysis/Shiny/Test/gaba_subset_forthpipe_03142023.RData")

# save dataset as a spare
gaba <- gaba.subset2

gaba@meta.data$all <- gaba@meta.data$seurat_clusters

inhib_clusters <- c("Pvalb","Sst","Sst Chodl","Meis/Ppp1r1b","Meis","Vip","Lamp5","removed")
inhib_hex <- (c('Pvalb'='#E66027','Sst'= '#F8991D',"Sst Chodl" = '#B0B235',"Meis/Ppp1r1b" = '#C03C82',
                "Meis" = '#C52126',"Vip" = '#A669AB',"Lamp5"='#DB808C',"removed"='#DB808C'))

levels(gaba@meta.data$all) <- inhib_clusters

for(c in 1:length(inhib_clusters)){
  gaba@meta.data$all[gaba@meta.data$all==(c-1)] <- inhib_clusters[c]
}

# creating a watermark text
watermark_text <- "Hope Lab©"  # Customize the watermark text

# Define UI for the Shiny app (side panel and main panel layout) ----
ui <- fluidPage(
  tags$head(
    tags$style(HTML(".selectize-dropdown { max-height: none; }"))
  ),
  # Theme for page -----
  theme=bs_theme(
    bootswatch = "minty"
  ),
  
  # App Title -----
  titlePanel("GABA Dataset"),
  
  # Sidebar layout with input and output definitions ---- 
  
  # selecting genes to get fed into plot ---- 
  
  sidebarLayout(
    sidebarPanel(
      style = "position:fixed;width:inherit;",
      width = 4,
      h3("Customize your plots!"),
      selectInput("Gr_UMAP",label="Choose to split UMAP by group or display all",
                  choices = c("Experiment Group" = "group",
                              "Rat" = "ratID",
                              "Sex" = "sex",
                              "All" = "all")),
      hr(),
      numericInput(
        "NoVars",
        "Number of genes to plot",
        value = 1,
        min = 1,
        max = 6
      ),
      uiOutput("VarsInput"),
      checkboxInput("show_ridge_clust", "Show cluster-specific Violin Plot"),
      conditionalPanel(
        condition = "input.show_ridge_clust == true",
        selectInput(
          "cluster",
          label = "Choose cluster to display",
          choices = c(inhib_clusters[1:6])
        )
      ),
      actionButton("submitBtn", "Submit", icon=icon("brain",stroke="black",class = "font-awesome",
                                                    style = "color: mediumvioletred;
                                                                   border-color: black"),
                   style= "color:black;
                                  background-color: powderblue")
      
      # adds a button to refresh plots by gene and group
    ),
    # Main panel for displaying outputs ----
    mainPanel(
      tabsetPanel(
        tabPanel("Single gene analysis",
                 plotOutput("umap_plot"),
                 fluidRow(
                   column(5,
                          plotOutput("umap_split"),
                          style='margin-bottom:25px;border:1px solid; padding: 5px;margin-left:55px;'
                   ),
                   column(5,
                          plotOutput("mainfeature"),
                          style='margin-bottom:25px;border:1px solid; padding: 5px;'
                   )
                 ),
                 hr(),
                 plotOutput("feature"),
                 hr(),
                 plotOutput("violin")
        ),
        
        # 2ND TAB ------
        tabPanel("Multiple gene analysis",
                 fluidRow(
                   column(6,
                          plotOutput("umap_plot2"),
                          style='margin-bottom:30px;border:1px solid; padding: 10px;'
                   ),
                   column(6,
                          plotOutput("mainfeature2"),
                          style='margin-bottom:30px;border:1px solid; padding: 10px;'
                   )
                 ),
                 hr(),
                 plotOutput("ridge2"),
                 hr(),
                 plotOutput("feature2")
        )
        
      )
    )
  )
)


# Define server logic for the Shiny app ---- 
server <- function(input, output, session){
  
  K <- reactive({
    input$NoVars
  })
  
  output$VarsInput <- renderUI({
    NoV <- K()
    Ge <- sapply(1:NoV, function(i) {
      paste0("gene", i)
    })
    Gr <- sapply(1:NoV, function(i) {
      paste0("group", i)
    })
    Gr_UMAP <- sapply(1:NoV, function(i) {
      paste0("groupu", i)
    })
    Gr_f <- sapply(1:NoV, function(i) {
      paste0("groupf", i)
    })
    
    output_list <- list()
    
    for (i in seq_along(1:NoV)) {
      output_list[[i]] <- fluidRow(
        column(
          4,
          textInput(
            Ge[i],
            label = "Choose gene",
            value = "",
            placeholder = "Gad1"
          )
        ),
        column(
          4,
          selectInput(
            Gr[i],
            label = "Split Violin Plot",
            choices = c(
              "Experiment" = "group",
              "Rat" = "ratID",
              "Sex" = "sex",
              "All" = "all"
            )
          ), class = "custom-dropdown"
        ),
        column(
          4,
          selectInput(
            Gr_f[i],
            label = "Split Feature Plot",
            choices = c(
              "Experiment" = "group",
              "Rat" = "ratID",
              "Sex" = "sex"
            )
          )
        )
      )
    }
    
    tagList(output_list)
  })
  
  # plots UMAP ----
  output$umap_plot <- renderPlot({
    umap <- DimPlot(gaba, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "all",
                    cols = inhib_hex, label.box=TRUE, label.color="white", repel=TRUE, label.size=5) +
      ggtitle(paste("Overall UMAP plot"))
    
    # Modify the plot to add borders and caption
    umap <- umap + 
      theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    # 
    # # Create the caption text
    # caption_text <- paste("Figure 1: UMAP Plot by ", input$Gr_UMAP, ".", sep="")
    # 
    # # Add the caption to the plot
    # umap <- umap + labs(caption = caption_text)
    # 
    
    # Display the modified plot
    print(umap)
    
  })
  
  # plots UMAP by group ----
  output$umap_split <- renderPlot({
    umap <- DimPlot(gaba, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = input$Gr_UMAP) +
      ggtitle(paste("UMAP plot by", input$Gr_UMAP))
    
    # Modify the plot to add borders and caption
    umap <- umap + 
      theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    # 
    # # Create the caption text
    # caption_text <- paste("Figure 1: UMAP Plot by ", input$Gr_UMAP, ".", sep="")
    # 
    # # Add the caption to the plot
    # umap <- umap + labs(caption = caption_text)
    # 
    
    # Display the modified plot
    print(umap)
    
  })
  
  
  
  # plots Violin ----
  output$violin <- renderPlot({
    plot_list <- list()
    
    for(i in 1:K()){
      gene <- input[[paste0("gene", i)]]
      group <- input[[paste0("group", i)]]
      
      validate(
        need(gene != "", "Please provide a gene.")
      )
      
      plot <- VlnPlot(gaba, features = str_to_title(gene), group.by = group) + 
        ggtitle(str_to_title(gene)) +
        theme(legend.position = "top", legend.justification = "center", plot.title = element_text(face = "bold.italic"))
      
      plot_list[[i]] <- plot
    }
    
    if (input$show_ridge_clust) {
      for(i in 1:K()){
        gene <- input[[paste0("gene", i)]]
        group <- input[[paste0("group", i)]]
        
        validate(
          need(gene != "", "Please provide a gene.")
        )
        
        plot <- VlnPlot(gaba, features = str_to_title(gene), group.by = group, idents = (which(inhib_clusters==input$cluster)-1)) + 
          ggtitle(str_to_title(gene)) +
          theme(legend.position = "top", legend.justification = "center", plot.title = element_text(face = "bold.italic"))
        
        plot_list[[i]] <- plot
      }
    }
    
    plot_output <- do.call(patchwork::wrap_plots, plot_list)
    plot_output
  })

  
  # plots Feature plot that is general ----
  output$mainfeature <- renderPlot({
    plot_list <- list()
    
    for(i in 1:K()){
      gene <- input[[paste0("gene", i)]]
      
      gene_i <- expr(italic(!!str_to_title(gene)))
      
      validate(
        need(gene != "", "Please provide a gene.")
      )
      
      plot <- FeaturePlot(gaba, features = str_to_title(gene)) +
        labs(title=expr(bold(paste("Feature Plot:",!!gene_i))))
      
      plot_list[[i]] <- plot
    }
    
    plot_output <- do.call(patchwork::wrap_plots, plot_list)
    print(plot_output)
  })
  
  # plots Feature plot split by group ----
  output$feature <- renderPlot({
    plot_list <- list()
    
    for(i in 1:K()){
      gene <- input[[paste0("gene", i)]]
      groupf <- input[[paste0("groupf", i)]]
      
      validate(
        need(gene != "", "Please provide a gene.")
      )
      
      plot <- FeaturePlot(gaba, features = str_to_title(gene), split.by = groupf) 
      
      plot_list[[i]] <- plot
    }
    
    plot_output <- do.call(patchwork::wrap_plots, plot_list)
    print(plot_output)
  })
}


# Run the Shiny app
shinyApp(ui = ui, server = server)
