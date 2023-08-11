# Figure 2 stats

# first calculate the percent per each rat

source("~/XPoSE/Scripts/Functions/calc_counts.R")

glut_group <- calc_counts(seur_obj = glut, fact1 = 'ratID',
                          fact2 = 'cluster_name',
                          fact3 = 'group', file_n = "group_glut.csv")

# combine by rat and group
glut_group$ratIDgroup <- paste0(glut_group$ratID, glut_group$group)

IDslist <- split(glut_group, glut_group[["ratIDgroup"]])

for (i in 1:length(IDslist)) {
        current_group <- IDslist[[i]]
        
        total_count <- sum(current_group[["count"]])
        current_group[["percent"]] <- current_group[["count"]] / total_count
        
        IDslist[[i]] <- current_group
}

# Combine the modified groups back into a single data frame
combined_df <- do.call(rbind, IDslist)
combined_df <- combined_df[complete.cases(combined_df[c("count", "percent")]), ]

# Z score -----------------------------------------------------------------

# Choose the reference group for calculating z-scores
reference_group <- "Homecage"  # Replace with the name of the reference group

clust_list <- split(combined_df, combined_df[["cluster_name"]])

anovas <- list()
anova_sum <- list()
posthoc_results <- list()

for (j in 1:length(clust_list)) {
        current_clust <- clust_list[[j]]
        
        # Subset the data for the reference group
        reference_subset <- subset(current_clust, group == reference_group)
        
        reference_mean <- mean(reference_subset[["percent"]])
        
        reference_sd <- sd(reference_subset[["percent"]])
        
        current_clust[["z_score"]] <- (current_clust[["percent"]] - reference_mean) / reference_sd
        
        #Run ANOVA
        current_aov <- clust_list[[j]]
        current_aov <- aov(z_score ~ group, data = current_clust)
        
        current_sum <- clust_list[[j]]
        current_sum <- summary(current_aov)
        
        current_posthoc <- clust_list[[j]]
        current_posthoc <- TukeyHSD(current_aov)
        
        clust_list[[j]] <- current_clust
        anovas[[j]] <- current_aov
        anova_sum[[j]] <- current_sum 
        #posthoc_results[[j]] <- current_posthoc
}

name <- names(clust_list)
names(anovas) <- name
names(anova_sum) <- name
#names(posthoc_results) <- name

combined_df2 <- do.call(rbind, clust_list)


