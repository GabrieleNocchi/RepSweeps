library(ggplot2)
library(viridis)

data <- readRDS("gab_picmin_results.rds")
data <- data$picmin_res

# Create categories based on FDR values
categories <- cut(data$picmin_fdr, breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1), labels = c("FDR < 0.1", "0.1 < FDR < 0.2", "0.2 < FDR < 0.3", "0.3 < FDR < 0.4", "0.4 < FDR < 0.5", "FDR > 0.5"))

# Convert to factor with desired order
categories <- factor(categories, levels = c("FDR < 0.1", "0.1 < FDR < 0.2", "0.2 < FDR < 0.3", "0.3 < FDR < 0.4", "0.4 < FDR < 0.5","FDR > 0.5"))

# Count occurrences in each category
count <- table(categories)

# Create a data frame for plotting
my_data <- data.frame(FDR = names(count), count)

# Plotting
p1 <- ggplot(my_data, aes(x = categories, y = count, fill = categories)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5, color = "black") +
  theme_bw() +
  ylab("OG count\n") +
  scale_fill_viridis_d(direction = 1, option = "mako") +
  ggtitle(paste("PicMin OGs =", format(round(as.numeric(length(data$p)), 1), big.mark = ","))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


p1b <- ggplot(my_data, aes(x = categories, y = count)) +
    geom_bar(stat = "identity", position = "stack", width = 0.5, color = "black", fill = "azure") +
    geom_text(aes(label = count, y = 5000), position = position_stack(vjust = 1.1), size = 4, color = "black") + # Specify y = 5000
    theme_bw() +
    ylab("OG count\n") +
    scale_fill_viridis_d(direction = 1, option = "mako") +
    ggtitle(paste("PicMin OGs =", format(round(as.numeric(length(data$p)), 1), big.mark = ","))) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    coord_flip() + theme_classic()

print(p1b)



# Assuming 'data' is your dataframe
library(ggplot2)

# Create a ggplot histogram
p2 <- ggplot(data, aes(x = picmin_fdr, fill = "PicMin FDR")) +
  geom_histogram(color="black", bins = 30, position = "identity", alpha = 0.7) +
  # Add labels and title
  labs(title = "PicMin FDR: 13,268 OGs",
       x = "PicMin FDR",
       y = "Frequency") +
  # Customize theme if needed
  theme_bw() +
  scale_fill_manual(values = "lightblue", name = "")



a <- readRDS("gab_picmin_results.rds")
a <- a$picmin_res
your_data <- a
# Assuming your data frame is named 'your_data'
library(ggplot2)

# Create the point plot with a more subdued color palette
p3 <- ggplot(your_data, aes(x = seq_along(Orthogroup), y = -log10(picmin_fdr), fill = as.factor(config_est))) +
  geom_point(shape = 21, size = 3, color = "black") + # Adding black point contours
  scale_fill_viridis_d() + # Using the viridis color palette
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red") + # Adding a red dashed line
  labs(x = "Orthogroups", y = "-log10(picmin_fdr)", fill = "Number of Species") + # Changing legend label
  theme_bw()+
  theme(
    legend.key.width = unit(0.2, "cm"), # Adjust the legend key width
    legend.key.height = unit(0.2, "cm") # Adjust the legend key height
  )



library(gridExtra)
grid.arrange(p1,p2,p3, ncol = 1)


gab_results <- readRDS("orthogroup_results.rds")

library(dplyr)

# Add a new column to your_data
your_data_2 <- your_data %>%
  mutate(count_ortho_DS_lt_0.1 = sapply(Orthogroup, function(og) sum(gab_results$Orthogroup == og & gab_results$ortho_DS < 0.1)))

p4 <- ggplot(your_data_2, aes(x = count_ortho_DS_lt_0.1, y = -log10(picmin_fdr))) +
  geom_point(shape = 21, size = 3, color = "black", fill = "grey") + # Adding black point contours
  scale_fill_viridis_d() + # Using the viridis color palette
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red") + # Adding a red dashed line
  labs(x = "Driving species (OmegaPlus emp-p < 0.1)", y = "-log10(picmin FDR)", fill = "Number of Species") + # Changing legend label
  theme_bw()
  
  
p5 <- ggplot(your_data_2, aes(x = count_ortho_DS_lt_0.1, y = -log10(picmin_fdr))) +
    geom_point(shape = 21, size = 3, color = "black", fill = "grey") + 
    scale_fill_viridis_d() + 
    geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red") + 
    ggtitle("OGs tested = 13,268") +
    labs(x = "Number of OP emp-p < 0.1", y = "-log10(picmin FDR)", fill = "Number of Species") + 
    theme_bw() +
    scale_x_continuous(breaks = seq(0, 8, by = 1)) +
    theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),axis.title = element_text(size = 10), title = element_text(size = 10))