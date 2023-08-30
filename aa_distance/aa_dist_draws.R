library(dplyr)

my_hits <- readRDS("gab_picmin_results.rds")

my_hits <- my_hits$picmin_res
to_crop <- my_hits$Orthogroup
to_crop <- as.data.frame(to_crop)
colnames(to_crop) <- "Orthogroup"
my_hits <- my_hits[my_hits$picmin_fdr < 0.5,]
my_hits <- my_hits$Orthogroup
my_hits <- as.data.frame(my_hits)
colnames(my_hits) <- "Orthogroup"


aa_dist <- read.table("formatted_distances.txt", h = F)
colnames(aa_dist) <- c("Orthogroup", "distance")

my_hits_dist <- my_hits %>%
  left_join(aa_dist, by = c("Orthogroup"))
my_hits_dist <-my_hits_dist[complete.cases(my_hits_dist), ]


aa_dist <- aa_dist %>%
  left_join(to_crop, by = c("Orthogroup"))
aa_dist <-aa_dist[complete.cases(aa_dist), ]

my_list <- replicate(10000,sample_n(aa_dist, length(my_hits_dist$Orthogroup)))

final_to_plot <- c()


for (i in 1:10000) {
    final_to_plot <- append(final_to_plot,mean(my_list[,i]$distance))
}


cat <- rep("All Distances", length(final_to_plot))

df <- cbind(final_to_plot,cat)



# Plot
df <- as.data.frame(df)

df$final_to_plot <- as.numeric(df$final_to_plot)

a <- mean(as.numeric(df$final_to_plot))
a1 <- quantile(df$final_to_plot,0.025)
a2 <- quantile(df$final_to_plot,0.975)


df <- cbind(a,"All Distances")
colnames(df) <- c("final_to_plot", "cat")


df <-as.data.frame(df)
df$final_to_plot <- as.numeric(df$final_to_plot)


library(ggplot2)
ggplot(df, aes(cat, (final_to_plot))) +        # ggplot2 plot with confidence intervals
    geom_point(color = "black", fill= "lightgoldenrod", shape = 19, size = 4) +
    geom_linerange(aes(x = "All Distances",ymin = a1, ymax = a2), size = 1.5) +
	 theme_classic() + coord_flip() +
    geom_point(aes(y=mean(my_hits_dist$distance), x = "All Distances"),colour="black", fill = "red", shape = 24, size = 3) +
    theme(axis.title.y=element_blank(),axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.title = element_text(hjust = 0.5, vjust = 1.5)) +
    ylab("Mean Orthogroup Amino Acid Distances") +  ggtitle("FDR < 0.5")
