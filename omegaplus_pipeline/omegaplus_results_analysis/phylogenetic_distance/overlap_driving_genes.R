library(ggplot2)
library(reshape2)
library(dplyr)
library(circlize)

z <- readRDS("gab_picmin_results.rds")
z <- z$picmin_res
z <- z[z$picmin_fdr<0.5,]
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

# Create the ggplot
ggplot(result, aes(x = phylo_dist, y = y_values)) +
  geom_point() +
  labs(x = "Phylogenetic distance", y = "Observed/Expected Hypergeometry Overlap of Driving Genes (emp-p < 0.1)") +
  ggtitle("Phylogenetic distance assessment") + theme_classic()
