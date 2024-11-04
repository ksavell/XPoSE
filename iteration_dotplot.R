# plot the enrichment DEG score

clusters <- unique(all$cluster_name)


all_data_lists <- list()

for (cl in clusters) {
  
  # Folder path is the same as cluster name in this case
  folder_id <- cl
  
  # Get a list of all CSV files in the folder
  csv_files <- list.files(path = cl, pattern = "results_.*_iteration_tally\\.csv", full.names = TRUE)
  
  # Initialize a list to store data frames for this specific folder
  data_list <- list()
  
  # Loop through each CSV file in the folder
  for (file in csv_files) {
    # Extract the unique number from the filename
    unique_number <- sub(".*results_(.*?)_iteration_tally\\.csv", "\\1", basename(file))
    
    # Create a unique name by combining the folder name and unique number
    data_name <- paste0(folder_id, "_", unique_number)
    
    # Read the CSV file into a data frame and store it in data_list
    data_list[[data_name]] <- read.csv(file)
  }

  # Store the data_list in all_data_lists with the cluster ID as the key
  all_data_lists[[folder_id]] <- data_list
}

# Check the structure to confirm
print(names(all_data_lists))

# now make the df for the plot --------------------------------------------

# Define clusters and numbers
clusters <- names(all_data_lists)
numbers <- c("0.035", "0.07", "0.105", "0.175", "0.245")

# Initialize data frames to store the results: one for averages, one for SD, and one for nods
average_sums_df <- data.frame(matrix(nrow = length(clusters), ncol = length(numbers)))
sd_sums_df <- data.frame(matrix(nrow = length(clusters), ncol = length(numbers)))

rownames(average_sums_df) <- clusters
colnames(average_sums_df) <- numbers
rownames(sd_sums_df) <- clusters
colnames(sd_sums_df) <- numbers

# Populate the data frames with mean, SD values, and nods data
for (cluster in clusters) {
  for (num in numbers) {
    # Construct the name of the DataFrame in all_data_lists
    df_name <- paste0(cluster, "_", num)
    
    # Check if the DataFrame exists and compute the average and SD if it does
    if (df_name %in% names(all_data_lists[[cluster]])) {
      df <- all_data_lists[[cluster]][[df_name]]
      numeric_df <- df[, sapply(df, function(col) is.numeric(col) || is.integer(col))]
      
      # Calculate column sums, then compute mean and SD
      col_sums <- colSums(numeric_df, na.rm = TRUE)
      average_sums_df[cluster, num] <- mean(col_sums)
      sd_sums_df[cluster, num] <- sd(col_sums)  # Calculate SD instead of SEM
    } else {
      # Assign NA if the DataFrame is missing
      average_sums_df[cluster, num] <- NA
      sd_sums_df[cluster, num] <- NA
    }
  }
}

# Check the resulting data frames
print(average_sums_df)
print(sd_sums_df)


# Add in the no downsample group ------------------------------------------

# read in the score column for each experience_tally file

# Define clusters and initialize an empty data frame to store the summed 'score' values
clusters <- names(all_data_lists)
scores_df <- data.frame(Cluster = character(), TotalScore = numeric(), stringsAsFactors = FALSE)

# Loop over each cluster and sum the 'score' column from experience_tally.csv
for (cl in clusters) {
  # Define the file path for experience_tally.csv in the current cluster folder
  experience_file <- file.path(cl, "experience_tally.csv")
  
  # Check if the experience_tally.csv file exists
  if (file.exists(experience_file)) {
    # Read the CSV file
    experience_data <- read.csv(experience_file)
    
    # Check if the 'score' column exists
    if ("score_column" %in% colnames(experience_data)) {
      # Calculate the sum of the 'score' column
      total_score <- sum(experience_data$score_column, na.rm = TRUE)
      
      # Add the total score to scores_df
      scores_df <- rbind(scores_df, data.frame(Cluster = cl, TotalScore = total_score))
    } else {
      warning(paste("No 'score' column found in", experience_file))
    }
  } else {
    warning(paste("experience_tally.csv not found for cluster:", cl))
  }
}

# View the combined scores_df
print(scores_df)

# add in columns that make it compatible with plotting
scores_df$Variable <- as.character(0.305)
scores_df$SD <- 0
scores_df <- scores_df[, c("Cluster", "Variable", "TotalScore", "SD")]
colnames(scores_df) <- c("Observation", "Variable", "Mean", "SD")

# dotplot time ------------------------------------------------------------

# Assuming average_sums_df and sd_sums_df are already created
library(ggplot2)
library(reshape2)

# Convert both average and SD data frames to long format
average_sums_df$Observation <- rownames(average_sums_df)
average_long <- melt(average_sums_df, id.vars = "Observation", variable.name = "Variable", value.name = "Mean")

sd_sums_df$Observation <- rownames(sd_sums_df)
sd_long <- melt(sd_sums_df, id.vars = "Observation", variable.name = "Variable", value.name = "SD")

# Merge the average and SD long data frames
plot_data <- merge(average_long, sd_long, by = c("Observation", "Variable"))

# Combine `scores_df` into `plot_data`
plot_data <- rbind(plot_data, scores_df)

# Define custom color mapping
custom_colors <- c(
  'CTL6' = '#2C8CB9',
  'PTL5' = '#0A5B8C',
  'ITL23' = '#41B75F',
  'ITL5' = '#5DBFC1', 
  'ITL6' = '#3A8F87',
  'NPL56' = '#3C9E64',
  'CTL6b' = '#6F499D',
  'Pvalb' = '#E66027',
  'Sst' = '#F8991D',
  'Meis2' = '#C52126',
  'Vip' = '#A669AB',
  'Lamp5' = '#DB808C',
  'SstChodl' = '#B0B235',
  'PvalbChand' = '#AD6C49'
)

# Plot with custom colors
ggplot(plot_data, aes(x = Variable, y = Mean, group = Observation, color = Observation)) +
  geom_line() +
  geom_point(data = subset(plot_data, !is.na(Mean)), size = 2) +  # Points for non-zero values
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0, alpha = 0.5, data = plot_data) +  # SD error bars
  labs(x = "Enrichment", y = "Mean ± SD") +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
