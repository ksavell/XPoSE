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
        
        # Calculate the total count within the current group
        total_count <- sum(current_group[["count"]])
        
        # Calculate the percentage for each cluster within the current group
        current_group[["percent"]] <- current_group[["count"]] / total_count
        
        # Replace the modified group back in the list
        IDslist[[i]] <- current_group
}

combined_df <- do.call(rbind, IDslist)

# Want to calculate z scores (homecage populations as the control)

# this I was testing from the manual proportion table I made in excel :O


anova_results <- aov(IT.L5.6 ~ group, data = glutT)

summary(anova_results)

posthoc_results <- TukeyHSD(anova_results)

print(posthoc_results)

