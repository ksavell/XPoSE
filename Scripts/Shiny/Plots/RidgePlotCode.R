#Ridge Plot Script (Both)
#By: Drake Thompson
#Date: Jun 2023
#Descript: Codes up the outputs for Ridge Plots for both GLUTamatergic
#          and GABA datasets

#Used libraries
library(shiny)
library(ggplot2)
library(Seurat)
library(umap)
library(dplyr)
library(patchwork)

#uses rel path to change directory
setwd("../data")
#getwd()

#stores loaded data as two separate vars
load("glut_subset_forthpipe_11102022.RData")
#glut <- load("glut_subset_forthpipe_11102022.RData")
load("gaba_subset_forthpipe_03142023.RData")

#gaba.subset2

# save dataset as a spare
gaba <- gaba.subset2
glut <- glut.subset2

##My Functions

#I am lazficient so I made a function so I don't have to run this twice
setDataUp <- function(dataSet){
    #Code I swiped from Padma that I turned into a function
    all.genes <- rownames(dataSet)
    DimPlot(dataSet, reduction = "umap", label = TRUE, pt.size = 0.5)
    
    cluster.ids <- names(dataSet$seurat_clusters)
    names(cluster.ids) <- levels(dataSet)
    gaba <- RenameIdents(dataSet, cluster.ids)
    DimPlot(dataSet, reduction = "umap", label = TRUE, pt.size = 0.5) # just an example of a plot
}


#Gets user input for settings
getWantedParams <- function(dataset){
    #Makes vector for feature names
    featVect <- dataset@assays[["RNA"]]@counts@Dimnames[[1]]
  
    #Gets feature(gene) name
    theFeat = readline(prompt = "What feature (gene symbol) would you like to plot?: ")
    
    #Converts to lowercase to ensure case insensitivity
    genesLower <- tolower(featVect)
    theFeat <- tolower(theFeat)
    
    #validates by checking if inputted name is within the lowercase feature list
    while (!(theFeat %in% genesLower)){
        #"error" message
        print("Invalid Input. Please return an existing gene symbol")
        #prompts user again
        theFeat = readline(prompt = "What feature would you like to plot?: ")
    }
    
    #Finds geneLower index that holds theFeat to convert back to original case
    index <- match(theFeat, genesLower)
    theFeat <- featVect[index]
    
    #Use to test if errors
    # testStr <- paste("This is the feat now: ", theFeat)
    # print(testStr)
        
    #creates prompt string
    promptStr <- "What would you like to group by?"
    opt1 <- "1. Group"
    opt2 <- "2. ratID"
    opt3 <- "3. sex"
    daPrompt <- paste(promptStr, opt1, opt2, opt3, sep = "\n")
    #prints prompt
    cat(daPrompt)
  
    #Gets user input
    theGroup = readline()
  
    vect <- c(1,2,3)  #Vector created to validate input.
    
    #validates input, loops if invalid
    while (!(theGroup %in% vect)) {
       #"error" message
       print("Invalid Input. Please return a 1, 2, or 3")
      
       #asks fo input again
       cat(daPrompt)
       theGroup = readline()
    }
    
    if (theGroup == 1){
      theGroup <- "Group"
    } else if (theGroup == 2){
      theGroup <- "ratID"
    } else {
      theGroup <- "sex"
    }
    
    #returns vars as a vector
    return(c(theFeat, theGroup))
}


#sets data up for both
#setDataUp(gaba)
#setDataUp(glut)

#Prints head and tail for testing purposes
print(head(glut.subset2@assays[["RNA"]]@counts@Dimnames[[1]]))
print(tail(glut.subset2@assays[["RNA"]]@counts@Dimnames[[1]]))

#Acquires params for object creation
#print("GABA Object Creation")
#gabaParams <- getWantedParams(gaba)
print("Glut Object Creation")
glutParams <- getWantedParams(glut)

#Finally makes plot, I included the ncol because I am unsure what happens without it
#RidgePlot(gaba, features = gabaParams[1], group.by = gabaParams[2], ncol = 2)
RidgePlot(glut, features = glutParams[1], group.by = glutParams[2])


 




#checks to see if loading was correct
#glut.subset2@assays$RNA
#gaba.subset2@assays$RNA



