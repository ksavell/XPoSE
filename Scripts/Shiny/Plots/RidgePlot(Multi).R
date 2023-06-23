#Ridge Plot Script Multi
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

#--------------------------------------------#
#YOU WILL LIKELY NEED TO CHNAGE THIS CODE

#uses rel path to change directory
setwd("../data")

#loads data
load("glut_subset_forthpipe_11102022.RData")
load("gaba_subset_forthpipe_03142023.RData")

#--------------------------------------------#

# save dataset into alias
gaba <- gaba.subset2
glut <- glut.subset2

#gets user's wanted dataset
promptStr <- "Which dataset would you like to use?"
opt1 <- "1. GABA"
opt2 <- "2. Glut"
daPrompt <- paste(promptStr, opt1, opt2, sep = "\n")
#prints prompt
cat(daPrompt)

#Gets user input
userData = readline()

vect <- c(1,2)  #Vector created to validate input.

#validates input, loops if invalid
while (!(userData %in% vect)) {
  #"error" message
  print("Invalid Input. Please return a 1 or 2.")
  
  #asks fo input again
  cat(daPrompt)
  userData = readline()
}

#stores dataset as general var
if (userData == 1){
  userData <- gaba
} else {
  userData <- glut
}

#Turns first letter to uppercase in most groups
userData@meta.data$sex[userData@meta.data$sex=="male"] <- "Male"
userData@meta.data$sex[userData@meta.data$sex=="female"] <- "Female"
userData@meta.data$group[userData@meta.data$group=="homecage"] <- "Homecage"
userData@meta.data$group[userData@meta.data$group=="negative"] <- "Negative"
userData@meta.data$group[userData@meta.data$group=="positive"] <- "Positive"


##My Functions


#validates input, asks user to reinput if not found
getFeature <- function(dataset, gene){
  #Makes vector for feature names
  featVect <- dataset@assays[["RNA"]]@counts@Dimnames[[1]]
  
  #Converts to lowercase to ensure case insensitivity
  genesLower <- tolower(featVect)
  theFeat <- tolower(gene)
  
  #validates by checking if inputted name is within the lowercase feature list
  while (!(theFeat %in% genesLower)){
    #"error" message
    errMessage <- paste("Invalid Input: `'", gene, "`' not found.", 
                        " Please return an existing gene symbol", "")
    print(errMessage)
    
    #prompts user again
    gene = readline(prompt = "What feature would you like to plot?: ")
    theFeat <- tolower(gene)
  }
  
  #Finds geneLower index that holds theFeat to convert back to original case
  index <- match(theFeat, genesLower)
  theFeat <- featVect[index]
  
  #returns correct feat
  return(theFeat)
}


#Gets user input for settings
getGroup <- function(dataset){
    
    #creates prompt string
    promptStr <- "What would you like to group by?"
    opt1 <- "1. Group"
    opt2 <- "2. ratID"
    opt3 <- "3. sex"
    opt4 <- "4. None"
    daPrompt <- paste(promptStr, opt1, opt2, opt3, opt4, sep = "\n")
    #prints prompt
    cat(daPrompt)
  
    #Gets user input
    theGroup = readline()
  
    vect <- c(1,2,3,4)  #Vector created to validate input.
    
    #validates input, loops if invalid
    while (!(theGroup %in% vect)) {
       #"error" message
       print("Invalid Input. Please return a 1, 2, 3 or 4.")
      
       #asks for input again
       cat(daPrompt)
       theGroup = readline()
    }
    
    #Converts group num to character
    if (theGroup == 1){
      theGroup <- "group"
    } else if (theGroup == 2){
      theGroup <- "ratID"
    } else if (theGroup == 3){
      theGroup <- "sex"
    }else {
      theGroup <- NULL
    }
    
    #returns corresponding group
    return(theGroup)
}


#sets data up for both
#setDataUp(gaba)
#setDataUp(glut)

#Prints head and tail for testing purposes
print(head(userData@assays[["RNA"]]@counts@Dimnames[[1]]))
print(tail(userData@assays[["RNA"]]@counts@Dimnames[[1]]))

#------------------------------------------------
#PADMA DELETE   PADMA DELETE    PADMA DELETE
inpVect <- c()
geneNum = readline("How many genes?: ")
geneNum <- as.numeric(geneNum)
for (i in 1:geneNum){
  gene <- readline("Input: ")
  inpVect <- append(inpVect, gene)
}

#It's stored as a vector, right?
#------------------------------------------------


#Acquires params for object creation
geneVect <- c()     #Empty vector to store (corrected) input

#loops for getFeature to clean user input
for (i in 1:length(inpVect)){
  #stores 'corrected' input
  gene <- getFeature(userData, inpVect[i])
  
  #checks for dupes, appends if not
  if (!(gene %in% geneVect)){
    geneVect <- append(geneVect, gene)
  }else {
    #doesn't append dupe
    print("Duplicate deleted")
  }
}

#swipes wanted group
theGroup <- getGroup(userData)

#Finally makes plot(s)
for (i in 1:length(geneVect)){
   #This print statement is mandatory (-_-)
   print(RidgePlot(userData, features = geneVect[i], group.by = theGroup))
}


 




#checks to see if loading was correct
#glut.subset2@assays$RNA
#gaba.subset2@assays$RNA



