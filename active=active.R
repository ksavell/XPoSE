# load packages needed
library(Seurat)
library(caret)
library(dplyr)
library(ggplot2)
library(devtools)
library(pROC)

# load Seurat object/data
load("~/Library/CloudStorage/Box-Box/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/DataForFigures/Robjects/glut.RData")

#glut_subset <- subset(glut, subset = group %in% c("Homecage", "Active", "Non-active"))

unique(glut@meta.data$group)

DimPlot(glut, group.by = "group")

glut$experience <- ifelse(glut@meta.data$group %in% c("Non-active", "Active"), "NC", "HC")

Idents(glut) <- glut$experience

roc <- FindMarkers(glut, ident.1 = "NC", ident.2 = "HC", test.use = "roc")

top_genes <- rownames(roc)[1:10]
top_genes_df <- as.data.frame(t(glut@assays$RNA@data[top_genes, ]))
top_genes_df$experience <- glut@meta.data$experience

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

full_features <- as.data.frame(t(GetAssayData(glut, layer = 'counts')[top_genes, ]))
predictions <- predict(model, full_features)

glut$predicted_group <- predictions

DimPlot(glut, group.by = "predicted_group")


# Augur package,  ---------------------------------------------------------
# why try this: is roc to find features appropriate? what is augur doing

obj <- glut

# extract the counts
obj_counts <- (obj[["RNA"]]$counts)

# now extract the metadata
obj_meta <- obj@meta.data

augur = calculate_auc(obj_counts, obj_meta, 
                      cell_type_col = "cluster_name", label_col = "experience",
                      rf_params = list(trees = 15, mtry = 100, min_n = NULL, importance = "accuracy")) # this line I got off a message board, using random forest

# plot the AUC values into umap
Idents(obj) <- "experience"
plot_umap(augur, obj, mode = "rank", palette = "Spectral", cell_type_col = "cluster_name") # this doesn't work.
