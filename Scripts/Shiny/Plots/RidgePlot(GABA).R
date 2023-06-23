#Ridge Plot Script (Glut)
#By: Drake Thompson
#Date: Jun 2023
#Descript: Codes up the outputs for Ridge Plots for GABA data

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
load("gaba_subset_forthpipe_03142023.RData")

#--------------------------------------------#


#shortens data name for ease of use
gaba <- gaba.subset2

#Makes vector for feature names
featVect <- gaba@assays[["RNA"]]@counts@Dimnames[[1]]

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
  theFeat <- tolower(theFeat)
}

#Finds geneLower index that holds theFeat to convert back to original case
index <- match(theFeat, genesLower)
theFeat <- featVect[index]

#creates prompt string for grouping name
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

#Converts num to group character
if (theGroup == 1){
  theGroup <- "Group"
} else if (theGroup == 2){
  theGroup <- "ratID"
} else {
  theGroup <- "sex"
}

#Finally makes plot
RidgePlot(gaba, features = theFeat, group.by = theGroup)









