library(ggplot2)
library(reshape2)
library(dplyr)
library(circlize)

z <- readRDS("gab_picmin_results.rds")
z <- z$picmin_res
z <- z[z$picmin_fdr < 0.5,]
z <- as.data.frame(cbind(z$Orthogroup, z$climate_var))
colnames(z) <- c("Orthogroup", "removeme")

df <- readRDS("orthogroup_results.rds")
df <- df %>% left_join(z, by = c("Orthogroup"))
df <- df[complete.cases(df),]
df <- df[, -c(2, 3, 4, 5, 6, 9, 10)]

df <- df[df$ortho_DS < 0.1,]

# Get unique Orthogroups in the dataframe
og_list <- unique(df$Orthogroup)

# Initialize a vector to store the mean of pairwise values
mean_pairwise_values <- vector("numeric", length = length(og_list))

# Initialize a vector to store the Orthogroup names
og_names <- vector("character", length = length(og_list))

# Iterate over each Orthogroup
for (i in 1:length(og_list)) {
  og <- og_list[i]

  # Subset the dataframe for the current Orthogroup
  og_df <- df[df$Orthogroup == og, ]

  species_list <- unique(og_df$species)
  comparison_file <- "all_phylo_dist.txt"

  # Read the comparison file
  comparison_df <- read.table(comparison_file, header = FALSE, stringsAsFactors = FALSE)

  # Initialize a list to store the pairwise values
  pairwise_values <- list()

  # Iterate over each row in the Orthogroup dataframe
  for (j in 1:(nrow(og_df) - 1)) {
    species1 <- og_df$species[j]

    # Iterate over the remaining species in the species_list
    for (k in (j + 1):nrow(og_df)) {
      species2 <- og_df$species[k]

      # Find the rows in the comparison file for the species pair
      row_indices <- which((comparison_df$V1 == species1 & comparison_df$V2 == species2) |
                            (comparison_df$V1 == species2 & comparison_df$V2 == species1))

      # If valid matches are found, store the pairwise comparison values
      if (length(row_indices) > 0) {
        pairwise_values[[length(pairwise_values) + 1]] <- comparison_df[row_indices, 3]
      }
    }
  }

  # Combine the pairwise values into a single vector
  pairwise_values <- unlist(pairwise_values)

  # Calculate the mean of pairwise values for the current Orthogroup
  mean_pairwise_values[i] <- mean(as.numeric(pairwise_values))

  # Store the Orthogroup name
  og_names[i] <- og
}

# Combine the Orthogroup names and mean pairwise values into a data frame
result_df <- data.frame(Orthogroup = og_names, Mean_Pairwise_Values = mean_pairwise_values)



# Load the ggplot2 library
library(ggplot2)

# Assuming 'result_df' is your data frame and 'Mean_Pairwise_Values' is the column containing the values you want to plot
ggplot(data = result_df, aes(x = Mean_Pairwise_Values)) +
  geom_histogram(bins = 100, fill = "gold", color = "black") +
  labs(title = "FDR < 0.5",
       x = "Mean pairwise phylogenetic distance of driving genes",
       y = "Frequency") + theme_classic()
