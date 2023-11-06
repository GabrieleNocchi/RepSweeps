library(ggplot2)
library(reshape2)
library(dplyr)
library(circlize)

z <- readRDS("gab_picmin_results.rds")
z <- z$picmin_res
z <- z[z$picmin_fdr<0.2,]
z <- as.data.frame(cbind(z$Orthogroup, z$climate_var))
colnames(z) <- c("Orthogroup", "removeme")


df <- readRDS("orthogroup_results.rds")
df <- df %>%   left_join(z, by = c("Orthogroup"))
df <- df[complete.cases(df),]
df <- df[, -c(2, 3, 4,5, 6, 9,10)]




library(dplyr)


# Step 2: Join the data with itself to get pairwise comparisons
pairwise_df <- df %>%
  inner_join(df, by = "Orthogroup") %>%
  filter(species.x != species.y) %>%
  group_by(species1 = species.x, species2 = species.y) %>%
  summarize(
    total_common_orthogroups = n_distinct(Orthogroup),
    species1_number_of_ortho_DS_less_than_0.1 = sum(ortho_DS.x < 0.1),
    species2_number_of_ortho_DS_less_than_0.1 = sum(ortho_DS.y < 0.1),
    shared_number_of_ortho_DS_less_than_0.1 = sum(ortho_DS.x < 0.1 & ortho_DS.y < 0.1)
  ) %>%
  ungroup()

# View the final dataframe
pairwise_df <- as.data.frame(pairwise_df)








# Create a new column with sorted species names
pairwise_df <- pairwise_df %>%
  mutate(sorted_species = ifelse(species1 < species2, paste(species1, species2, sep = "_"), paste(species2, species1, sep = "_")))

# Remove duplicate rows based on the sorted_species column

pairwise_df <- pairwise_df %>%
  distinct(sorted_species, .keep_all = TRUE) %>%
  select(-sorted_species)


all_dist <- read.table("all_phylo_dist.txt", h = T)


# Perform the inner join based on species1 and species2 columns
result <- inner_join(pairwise_df, all_dist, by = c("species1", "species2"))

# Plot


# Load the required libraries
library(ggplot2)

# Assuming you already have 'result' dataframe from the previous step

# Calculate the values for the y-axis
y_values <- result$shared_number_of_ortho_DS_less_than_0.1 / ((result$species1_number_of_ortho_DS_less_than_0.1 * result$species2_number_of_ortho_DS_less_than_0.1) / result$total_common_orthogroups)

cor_1 <- cor.test(result$phylo_dist, y_values)

# Create the ggplot
p1 <- ggplot(result, aes(x = phylo_dist, y = y_values)) +
    geom_point(color = "black", cex  = 0.5) +  # Changed point color to red
    geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue", alpha = 0.5) +  # Added correlation line with shading
    labs(x = "Phylogenetic distance", y = "Observed/Expected", cex = 2) +
    ggtitle(paste("FDR < 0.2 (",round(cor_1$estimate,2),", p-value=", round(cor_1$p.value,2),")"))  +
    scale_color_manual(values = "black")  + theme_bw() + # Changed point color to red
    theme(
      axis.title.x = element_text(size =8),  # Adjust x-axis label size
      axis.title.y = element_text(size = 8),  # Adjust y-axis label size
      plot.title = element_text(size = 8)   # Adjust plot title size
    )







    library(ggplot2)
    library(reshape2)
    library(dplyr)
    library(circlize)

    z <- readRDS("gab_picmin_results.rds")
    z <- z$picmin_res
    z <- z[z$picmin_fdr<0.1,]
    z <- as.data.frame(cbind(z$Orthogroup, z$climate_var))
    colnames(z) <- c("Orthogroup", "removeme")


    df <- readRDS("orthogroup_results.rds")
    df <- df %>%   left_join(z, by = c("Orthogroup"))
    df <- df[complete.cases(df),]
    df <- df[, -c(2, 3, 4,5, 6, 9,10)]




    library(dplyr)


    # Step 2: Join the data with itself to get pairwise comparisons
    pairwise_df <- df %>%
      inner_join(df, by = "Orthogroup") %>%
      filter(species.x != species.y) %>%
      group_by(species1 = species.x, species2 = species.y) %>%
      summarize(
        total_common_orthogroups = n_distinct(Orthogroup),
        species1_number_of_ortho_DS_less_than_0.1 = sum(ortho_DS.x < 0.1),
        species2_number_of_ortho_DS_less_than_0.1 = sum(ortho_DS.y < 0.1),
        shared_number_of_ortho_DS_less_than_0.1 = sum(ortho_DS.x < 0.1 & ortho_DS.y < 0.1)
      ) %>%
      ungroup()

    # View the final dataframe
    pairwise_df <- as.data.frame(pairwise_df)








    # Create a new column with sorted species names
    pairwise_df <- pairwise_df %>%
      mutate(sorted_species = ifelse(species1 < species2, paste(species1, species2, sep = "_"), paste(species2, species1, sep = "_")))

    # Remove duplicate rows based on the sorted_species column

    pairwise_df <- pairwise_df %>%
      distinct(sorted_species, .keep_all = TRUE) %>%
      select(-sorted_species)


    all_dist <- read.table("all_phylo_dist.txt", h = T)


    # Perform the inner join based on species1 and species2 columns
    resulted <- inner_join(pairwise_df, all_dist, by = c("species1", "species2"))

    # Plot


    # Load the required libraries
    library(ggplot2)

    # Assuming you already have 'result' dataframe from the previous step

    # Calculate the values for the y-axis
    y_valuess <- resulted$shared_number_of_ortho_DS_less_than_0.1 / ((resulted$species1_number_of_ortho_DS_less_than_0.1 * resulted$species2_number_of_ortho_DS_less_than_0.1) / resulted$total_common_orthogroups)
	
	cor_2 <- cor.test(result$phylo_dist, y_valuess)

    # Create the ggplot
    p2 <- ggplot(resulted, aes(x = phylo_dist, y = y_valuess)) +
        geom_point(color = "black", cex = 0.5) +  # Changed point color to red
        geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue", alpha = 0.5) +  # Added correlation line with shading
        labs(x = "Phylogenetic distance", y = "Observed/Expected") +
        ggtitle(paste("FDR < 0.1 (",round(cor_2$estimate,2),", p-value=", round(cor_2$p.value,2),")")) +
        scale_color_manual(values = "black")  + theme_bw() + # Changed point color to red
        theme(
          axis.title.x = element_text(size =8),  # Adjust x-axis label size
          axis.title.y = element_text(size = 8),  # Adjust y-axis label size
          plot.title = element_text(size = 8)   # Adjust plot title size
        )


setwd(dir = "./mean_OG_phylogenetic_distance/")



library(ggplot2)
library(reshape2)
library(dplyr)
library(circlize)

z <- readRDS("gab_picmin_results.rds")
z <- z$picmin_res
z <- z[z$picmin_fdr < 0.2,]
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
p3 <- ggplot(data = result_df, aes(x = Mean_Pairwise_Values)) +
  geom_histogram(bins = 100, fill = "gold", color = "black") +
  labs(title = "FDR < 0.2",
       x = "Mean phylogenetic distance of driving genes",
       y = "Frequency") + theme_bw() +
       theme(
         axis.title.x = element_text(size =8),  # Adjust x-axis label size
         axis.title.y = element_text(size = 8),  # Adjust y-axis label size
         plot.title = element_text(size = 8)   # Adjust plot title size
       )








       library(ggplot2)
       library(reshape2)
       library(dplyr)
       library(circlize)

       z <- readRDS("gab_picmin_results.rds")
       z <- z$picmin_res
       z <- z[z$picmin_fdr < 0.1,]
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
       resulted_df <- data.frame(Orthogroup = og_names, Mean_Pairwise_Values = mean_pairwise_values)



       # Load the ggplot2 library
       library(ggplot2)

       # Assuming 'result_df' is your data frame and 'Mean_Pairwise_Values' is the column containing the values you want to plot
      p4 <-  ggplot(data = resulted_df, aes(x = Mean_Pairwise_Values)) +
         geom_histogram(bins = 100, fill = "gold", color = "black") +
         labs(title = "FDR < 0.1",
              x = "Mean phylogenetic distance of driving genes",
              y = "Frequency") + theme_bw() +
              theme(
                axis.title.x = element_text(size =8),  # Adjust x-axis label size
                axis.title.y = element_text(size = 8),  # Adjust y-axis label size
                plot.title = element_text(size = 8)   # Adjust plot title size
              )



library(gridExtra)
grid.arrange(p1,p2,p3,p4,ncol =2)
