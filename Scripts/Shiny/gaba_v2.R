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

setwd("~/Desktop/Analysis/Shiny/Test")

load("~/Desktop/Analysis/Shiny/Test/gaba_subset_forthpipe_03142023.RData")

# save dataset as a spare
gaba <- gaba.subset2

gaba@meta.data$all <- gaba@meta.data$seurat_clusters

# creating a watermark text
watermark_text <- "Hope Lab©"  # Customize the watermark text

# Define UI for the Shiny app (side panel and main panel layout)
ui <- fluidPage(
  
  # Theme for page -----
  theme=bs_theme(
    bootswatch = "minty"
  ),
  
  # App Title -----
  titlePanel("GABA Dataset"),
  
  # Sidebar layout with input and output definitions ---- 
 # sidebarLayout(
    
    # selecting genes to get fed into plot ---- 
   
               sidebarLayout(
               # fluidRow(
               #   column(
               #     9,
                   sidebarPanel(
                     h3("Customize your plots!"),
                     numericInput(
                       "NoVars",
                       "Number of genes to plot",
                       value = 1,
                       min = 1,
                       max = 6
                     ),
                     uiOutput("VarsInput")
                     
                     # adds a button to refresh plots by gene and group
                     #submitButton("Create plots",icon=fa("brain",fill="pink",stroke="black"))
                   ),
               #    column(3,
                          # Main panel for displaying outputs ----
                          mainPanel(
                            tabsetPanel(
                              tabPanel("GABA Dataset",
                            #plotOutput("plot"),
                            plotOutput("umap_plot"),
                            plotOutput("ridge"),
                            plotOutput("feature")
                            
                            # fluidRow(
                            #   column(6,
                            #          plotOutput("ridge")),
                            #   column(6,
                            #          plotOutput("feature"))
                            # ))
               # )
                   
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
    NoV = K()
    Ge = sapply(1:NoV, function(i){paste0("gene",i)})
    Gr = sapply(1:NoV, function(i){paste0("group",i)})
    Gr_UMAP = sapply(1:NoV, function(i){paste0("groupu",i)})
    
    output = tagList()
    
    for(i in seq_along(1:NoV)){
      output[[i]] = tagList()
      output[[i]][[1]] = textInput(Ge[i], label = "Choose gene", value = "", placeholder = "Gad1")
      output[[i]][[2]] =  selectInput(Gr_UMAP[i],label="Choose to split UMAP by group or display all",
                                      choices = c("Experiment Group" = "group",
                                                  "Rat" = "ratID",
                                                  "Sex" = "sex",
                                                  "All" = "all"))
      output[[i]][[3]] =  selectInput(Gr[i],label="Choose to split other plots by group or display all",
                                      choices = c("Experiment Group" = "group",
                                                  "Rat" = "ratID",
                                                  "Sex" = "sex",
                                                  "All" = "all"))
    }
    
    output
  })
  
  # # Sanitizes the error message for the UI
  # gene <- reactive({
  #   validate(
  #     need(input$gene != "", "Please provide a gene.")
  #   )
  # })
  
  # plots UMAP ----
  output$umap_plot <- renderPlot({
    #  umap <- DimPlot(gaba, reduction = "umap", label = TRUE, pt.size = 0.5) -- his line works fine
    # print(umap)
    plot_list <- list()
    # gaba.features <- c("Gad1","Kcnc2","Sst","Chodl","Ppp1r1b","Meis2","Vip","Lamp5")
    for(i in 1:K()){
      groupu <- input[[paste0("groupu", i)]]
      
      plot <- DimPlot(gaba, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = groupu) +
        ggtitle(NULL)
      
      plot_list[[i]] <- plot
    }
    
    plot_output <- do.call(patchwork::wrap_plots, plot_list) 
    print(plot_output)
  })
  
  
  
  # plots Ridge ----

  output$ridge <- renderPlot({
    plot_list <- list()
    
    for(i in 1:K()){
      gene <- input[[paste0("gene", i)]]
      group <- input[[paste0("group", i)]]
      
      validate(
        need(gene != "", "Please provide a gene.")
      )
      
      plot <- RidgePlot(gaba, features = str_to_title(gene), group.by = group) + 
        theme(legend.position = "top", legend.justification = "center")
      
      plot_list[[i]] <- plot
    }
    
    plot_output <- do.call(patchwork::wrap_plots, plot_list)
    print(plot_output)
    
  })
  
  # plots Feature plot ----
  output$feature <- renderPlot({
    plot_list <- list()
    
    for(i in 1:K()){
      gene <- input[[paste0("gene", i)]]
      group <- input[[paste0("group", i)]]
      
      validate(
        need(gene != "", "Please provide a gene.")
      )
      
      plot <- FeaturePlot(gaba, features = str_to_title(gene), split.by = group) 
      
      plot_list[[i]] <- plot
    }
    
    plot_output <- do.call(patchwork::wrap_plots, plot_list)
    print(plot_output)
    # 
    # feature <- FeaturePlot(gaba, features = str_to_title(C))
    # vals$feature <- feature
    # print(feature)
  })
  
  # Downloadable plot of selected gene and group ----
  # 
  # output$downloadUMAP <-  downloadHandler( # downloads UMAP
  #   filename = function() { 
  #     paste("plot_", "UMAP",".png",sep="")
  #   },
  #   
  #   content = function(file) {
  #     png(file, width=500, height=500)
  #     print(vals$umap)
  #     dev.off()
  #     
  #     # Open the generated plot
  #     img <- image_read(file)
  #     
  #     # Add the watermark to the plot
  #     img_with_watermark <- image_annotate(img, watermark_text, location = "+300+150", color = "grey", size = 20)
  #     
  #     # Save the modified plot with the watermark
  #     image_write(img_with_watermark, path = file)
  #   }
  # )
  # 
  # output$downloadRidge <- downloadHandler( # downloads Ridge Plot
  #   filename = function(){
  #     paste("plot_",str_to_title(input$gene), "RidgePlot_by",input$group,".png",sep="")
  #   },
  #   content = function(file) {
  #     png(file,width=1000)
  #     print(vals$ridge)
  #     dev.off()
  #   }
  # )
  # 
  # output$downloadFeature <- downloadHandler( # downloads Feature Plot
  #   filename = function(){
  #     paste("plot_",str_to_title(input$gene), "FeaturePlot",".png",sep="")
  #   },
  #   content = function(file) {
  #     png(file,width=1000,height=1000)
  #     print(vals$feature)
  #     dev.off()
  #   }
  # )
  # 
  
  
  # output$downloadData <-  downloadHandler(
  #   filename = function() {
  #     paste("plot_", input$downloadplot, ".pdf",sep="")
  #   },
  #   content = function(file){
  #     pdf(file,width=15,height=8)
  #     print(vals$umap)
  #     dev.off()
  #   })
}


# Run the Shiny app
shinyApp(ui = ui, server = server)
