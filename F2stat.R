# Figure 2 stats

# first calculate the percent per each rat

source("~/XPoSE/Scripts/Functions/calc_counts.R")

glut_group <- calc_counts(seur_obj = glut, fact1 = 'ratID',
                          fact2 = 'cluster_name',
                          fact3 = 'group', file_n = "group_glut.csv")

# combine by rat and group
glut_group$ratIDgroup <- paste0(glut_group$ratID, glut_group$group)

# split by rat and find percent for each cluster per population
IDslist <- split(glut_group, glut_group[["ratIDgroup"]])




calc_proportion <- function(countdf, split_fact = "ratID") {
        IDslist <- split(countdf, countdf[[split_fact]])
        
        IDslist <- sapply(IDslist, function(x) {
                x$prop <- prop.table(table(x$count))
        })
        
        return(IDslist)
}

calc_proportion <- function(countdf, split_fact = "ratIDgroup") {
        # Split the data frame by 'split_fact' (ratIDgroup)
        IDslist <- split(countdf, countdf[[split_fact]])
        
        calc_group_proportions <- function(group_df) {
                proportions <- prop.table(table(group_df$cluster_name, group_df$count)) # this might be where it's going wrong
                prop_df <- as.data.frame(proportions)
                prop_df <- cbind(cluster_name = rownames(prop_df), prop_df)
                return(prop_df)
        }
        
        # Calculate proportions for each group
        prop_list <- lapply(IDslist, calc_group_proportions)
        
        return(prop_list)
}

test <- calc_proportion(glut_group, split_fact = "ratIDgroup")


# Want to calculate z scores (homecage populations as the control)

# this I was testing from the manual proportion table I made in excel :O


anova_results <- aov(IT.L5.6 ~ group, data = glutT)

summary(anova_results)

posthoc_results <- TukeyHSD(anova_results)

print(posthoc_results)

