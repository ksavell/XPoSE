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

setwd("~/Desktop/Analysis/Shiny/Test")

load("/Volumes/CM-Data2/Robjects/gaba_subset_forthpipe_11102022.RData")

# save dataset as a spare
gaba <- gaba.subset2

# Define UI for the Shiny app (side panel and main panel layout)
ui <- fluidPage(
  
  # Theme for page -----
  theme=bs_theme(
    bootswatch = "morph"
  ),
  
  # App Title -----
  titlePanel("GABA Datasets"),
  
  # Sidebar layout with input and output definitions ---- 
  sidebarLayout(
    
    # Sidebar panel for inputs ----- 
    sidebarPanel(
      h3("Customize your plots!"),
      # adds user input box for gene(s) of choice
      textInput("gene",
                label="Choose a gene",
                value="",
                placeholder="Gad1"),
      
      # adds user input box for group of choice or by all
      selectInput("group",
                  label="Choose to split by group or display all",
                  choices=c("Experiment Group" = "group",
                            "Rat" = "ratID", 
                            "Sex" = "sex",
                            "All")),
      
      # adds a button to refresh plots by gene and group
      submitButton("Create plots",icon=fa("brain",fill="pink",stroke="black")),
      
      br(), #adds a line break
      # 
      # # adds a download format button
      # radioButtons("format",
      #              label="Download format",
      #              choices=c("PDF"=".pdf",
      #                        "PNG"=".png")),
      # 
      # # adds a download button and chooses which plot to download
      # radioButtons("downloadplot",label="Choose which plot to download",
      #             choices = c("UMAP",
      #                         "Ridge Plot")),
      
      downloadButton("downloadUMAP", "Download UMAP as PDF", icon = shiny::icon("download")),
      br(),
      br(),
      downloadButton("downloadRidge", "Download Ridge Plot as PDF", icon = shiny::icon("download"))
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      #plotOutput("plot"),
      plotOutput("umap_plot"),
      plotOutput("ridge")
    )
    
  )
)

# Define server logic for the Shiny app ---- 
server <- function(input, output, session){
  
  vals <- reactiveValues()
  
  # Sanitizes the error message for the UI
  gene <- reactive({
    validate(
      need(input$gene != "", "Please provide a gene.")
    )
  })
  
  # plots UMAP ----
  output$umap_plot <- renderPlot({
    umap <- DimPlot(gaba, reduction = "umap", label = TRUE, pt.size = 0.5)
    vals$umap <- umap
    print(umap)
  })
  
  
  
  # plots Ridge ----
  # manipulates stored variable if user chooses not to group
  g <- reactive({
    if(input$group=="All"){
      NULL
    }
    else{
      input$group
    }
  })
  
  output$ridge <- renderPlot({
    head(gene())
    ridge_final <- RidgePlot(gaba, features = str_to_title(input$gene), group.by = g())
    vals$ridge <- ridge_final
    print(ridge_final)
  })
  
  # plot <- reactive({
  #   if(input$downloadplot=="UMAP"){
  #     output$umap_plot
  #   }
  #   else if(input$downloadplot=="Ridge plot"){
  #     output$ridge
  #   }
  # })
  # 
  
  # Downloadable plot of selected gene and group ----
  
  output$downloadUMAP <-  downloadHandler(
    filename = function() {
      paste("plot_", "UMAP",".pdf",sep="")
    },
    
    content = function(file) {
      pdf(file,width=15,height=8)
      print(vals$umap)
      dev.off()
    }
  )
  
  output$downloadRidge <- downloadHandler(
    filename = function(){
      paste("plot_",str_to_title(input$gene), "RidgePlot_by",input$group,".pdf",sep="")
    },
    content = function(file) {
      pdf(file,width=15,height=8)
      print(vals$ridge)
      dev.off()
    }
  )
  
  
  
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
