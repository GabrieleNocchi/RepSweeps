a <- readRDS("gab_picmin_results.rds")
a <- a$picmin_res
your_data <- a
# Assuming your data frame is named 'your_data'
library(ggplot2)

# Create the point plot with a more subdued color palette
ggplot(your_data, aes(x = seq_along(Orthogroup), y = -log10(picmin_fdr), fill = as.factor(config_est))) +
  geom_point(shape = 21, size = 3, color = "black") + # Adding black point contours
  scale_fill_viridis_d() + # Using the viridis color palette
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red") + # Adding a red dashed line
  labs(x = "Orthogroups", y = "-log10(picmin_fdr)", fill = "Number of species driving repeatability") + # Changing legend label
  theme_bw()
