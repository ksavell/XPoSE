# load packages needed
library(Seurat)
library(caret)
library(dplyr)
library(ggplot2)
library(devtools)
library(pROC)

# load Seurat object/data
load("/Users/holmesar/Library/CloudStorage/Box-Box/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/DataForFigures/Robjects/glut.RData")

glut_subset <- subset(glut, subset = group %in% c("Homecage", "Active", "Non-active"))

unique(glut_subset@meta.data$group)

DimPlot(glut_subset, group.by = "group")

glut_subset$experience <- ifelse(glut_subset@meta.data$group %in% c("Non-active", "Homecage"), "Control", "Active")

Idents(glut_subset) <- glut_subset$experience

roc <- FindMarkers(glut_subset, ident.1 = "Active", ident.2 = "Control", test.use = "roc")

top_genes <- rownames(roc)[1:10]
top_genes_df <- as.data.frame(t(glut_subset@assays$RNA@data[top_genes, ]))
top_genes_df$experience <- glut_subset@meta.data$experience

set.seed(150)
trainIndex <- createDataPartition(top_genes_df$experience, p = 0.8, list = FALSE)
train_data <- top_genes_df[trainIndex, ]
test_data <- top_genes_df[-trainIndex, ]

train_features <- train_data %>% select(-experience)
train_labels <- train_data$experience
test_features <- test_data %>% select(-experience)
test_labels <- test_data$experience

model <- train(train_features, train_labels, method = "rf",
               trControl = trainControl(method = "cv", number = 5))

predictions <- predict(model, newdata = test_features)

predictions <- as.factor(predictions)
test_labels <- as.factor(test_labels)

conf_matrix <- confusionMatrix(predictions, test_labels)
print(conf_matrix)

accuracy <- sum(predictions == test_labels) / length(test_labels)
print(paste("Accuracy:", accuracy))

full_features <- as.data.frame(t(GetAssayData(glut_subset, slot = 'counts')[top_genes, ]))
predictions <- predict(model, full_features)

glut_subset$predicted_group <- predictions

DimPlot(glut_subset, group.by = "predicted_group")
