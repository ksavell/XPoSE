#Makes a volcano plot for the data given
#Authors: Katherine Savell, Padma Saravanan, Drake Thompson

# Info --------------------------------------------------------------------

#This script currently uses a tibble made from other tibbles w/ the GABA and glut
#datasets. It may have to be generated prior to use or we can turn it to a csv at
#a later date.


# Packages ----------------------------------------------------------------
library(shiny)
library(plotly)
library(dplyr)
library(tidyverse)


# Data Treating -----------------------------------------------------------
#makes data var
df <- coexp     #You will need to generate this from the upset plot code

#Labels different parts
df$diffexpressed0[df$log2FoldChange_0 > 0 & df$padj_0 < 0.05] <- "UP"
df$diffexpressed0[df$log2FoldChange_0 < 0 & df$padj_0 < 0.05] <- "DOWN"
df$delabel<- NA

#makes new columns so hover display works
df$`Log2 (fold change)` <- as.numeric(format(round(df$log2FoldChange_0, 3))) 
df$`Adjusted p-value`<- -log10(df$padj_0)


# Plot Creation -----------------------------------------------------------

#makes plot
volc_plot <- ggplot(data=(df),
          aes(x=`Log2 (fold change)`,
              y=`Adjusted p-value`,
              col=diffexpressed0,
              label=delabel))+
  geom_point(size= 3.0)+
  theme_classic()+
  geom_text()+
  #Cool colors that have some degree of contrast
  scale_color_manual(values=c('darkorchid','darkorange', 'transparent'))+
  theme(text=element_text(size=16))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0, size = 1.5)+
  #labels axes
  xlab("Log2 (fold change)")+
  ylab("Adjusted Pvalue") +
  #Range- should be changed with different data
  xlim(c(-2,2))

#Makes plot interactable
volc_plot <- ggplotly(volc_plot, tooltip = c("x", "y"))


# Shiny Stuff -------------------------------------------------------------
shinyApp(
  ui <- shinyUI(fluidPage(
    sidebarLayout(sidebarPanel( h4("Test Plot")),
                  mainPanel(plotlyOutput("plot1"))
    )
  )),
  
  server <- shinyServer(
    function(input, output) {
      output$plot1 <- renderPlotly({
        volc_plot 
      })
    }
  ))

shinyApp(ui, server)
