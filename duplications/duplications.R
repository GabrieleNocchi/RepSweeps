# Read the Tau scores into my_data
# map is what I use to link TAIR_gene to Orthogroup_ID
my_hits <- readRDS("gab_picmin_results.rds")
my_hits <- my_hits$picmin_res
my_hits <- my_hits[my_hits$picmin_fdr < 0.5,]
my_hits <- my_hits$Orthogroup
my_hits <- as.data.frame(my_hits)
colnames(my_hits) <- "Orthogroup"


to_crop <- readRDS("gab_picmin_results.rds")
to_crop <- to_crop$picmin_res


map <- read.table("all_duplications.txt", h = T)

library(dplyr)
all <- to_crop %>%
  left_join(map, by = c("Orthogroup"))
  all<-all[complete.cases(all), ]


my_hits_z <- my_hits %>%
  left_join(all, by = c("Orthogroup"))
my_hits_z1 <-my_hits_z[complete.cases(my_hits_z), ]


### DRAWS

my_list <- replicate(10000,sample_n(all, length(my_hits_z1$Orthogroup)))

final_to_plot <- c()


for (i in 1:10000) {
    final_to_plot <- append(final_to_plot,mean(my_list[,i]$occurrences))
}


cat <- rep("All Duplications", length(final_to_plot))

df <- cbind(final_to_plot,cat)

# Plot
df <- as.data.frame(df)

df$final_to_plot <- as.numeric(df$final_to_plot)

a <- mean(as.numeric(df$final_to_plot))
a1 <- quantile(df$final_to_plot,0.025)
a2 <- quantile(df$final_to_plot,0.975)


df <- cbind(a,"All Duplications")
colnames(df) <- c("final_to_plot", "cat")


df <-as.data.frame(df)
df$final_to_plot <- as.numeric(df$final_to_plot)


library(ggplot2)
ggplot(df, aes(cat, (final_to_plot))) +        # ggplot2 plot with confidence intervals
    geom_point(color = "black", fill= "lightgoldenrod", shape = 19, size = 4) +
    geom_linerange(aes(x = "All Duplications",ymin = a1, ymax = a2), size = 1.5) +
	 theme_classic() + coord_flip() +
    geom_point(aes(y=mean(my_hits_z1$occurrences), x = "All Duplications"),colour="black", fill = "red", shape = 24, size = 3) +
    theme(axis.title.y=element_blank(),axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.title = element_text(hjust = 0.5, vjust = 1.5)) +
    ylab("Mean Orthogroup Duplications") +  ggtitle("FDR < 0.5")
