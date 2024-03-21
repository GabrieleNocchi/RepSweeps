library(dplyr)
library(ggplot2)
library(gridExtra)

results <- readRDS("orthogroup_results.rds")
results$difference <- sapply(strsplit(sub(".*:", "", results$gene), "-"), function(x) as.numeric(x[2]) - as.numeric(x[1]))

unique_species <- unique(results$species)

# Initialize an empty list to store plots
plots_list <- list()
all_corr <- data.frame()

for (specie in unique_species) {
  species_data <- filter(results, species == specie)

  # Calculate correlation
  cor_test_result <- cor.test(species_data$difference, -log10(species_data$mean_emp_p))

  # Extract correlation coefficient and p-value
  cor_coef <- cor_test_result$estimate
  p_value <- cor_test_result$p.value
  my_row <- cbind(specie, mean(species_data$difference), cor_coef)
  all_corr <- rbind(all_corr,my_row)

  # Create the scatter plot with title
  current_plot <- ggplot(species_data, aes(x = difference, y = -log10(mean_emp_p))) +
    geom_point(size = 0.5, shape = 21, fill = "yellow") +
    labs(title = paste(specie," Mean gene length:",round(mean(species_data$difference),2), "\n", round(cor_coef, 2), "P-value:", format.pval(p_value)),
         x = "Gene Length",
         y = "-log(10) emp-p")+
  theme(
    plot.title = element_text(size = 6, face = "bold"),  # Adjust title size and style
    axis.title.x = element_text(size = 5),  # Adjust x-axis label size
    axis.title.y = element_text(size = 5),  # Adjust y-axis label size
    axis.text.x = element_text(size = 5),   # Adjust x-axis tick label size
    axis.text.y = element_text(size = 5)    # Adjust y-axis tick label size
  )

  # Append the current plot to the list
  plots_list[[length(plots_list) + 1]] <- current_plot
}
svg("length.svg", height = 12,width = 8)
# Arrange the plots using grid.arrange
grid.arrange(grobs = plots_list, ncol = 4)  # Adjust ncol as needed
dev.off()



corr_2 <- cor.test(as.numeric(all_corr$V2),as.numeric(all_corr$cor_coef))
ggplot(all_corr, aes(x = round(as.numeric(V2),2), y = as.numeric(cor_coef))) +
  geom_point(size = 5, shape = 21, fill = "yellow")+
  labs(title = paste("Correlation", "\n", round(corr_2$estimate, 2), "P-value:", format.pval(corr_2$p.value)),
       x = "Species Average Gene Length",
       y = "correlation between gene size and -log(10emp-p)")
