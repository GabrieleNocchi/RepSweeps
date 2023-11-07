library(dplyr)
library(ggplot2)

my_hits <- readRDS("gab_picmin_results.rds")

my_hits <- my_hits$picmin_res
to_crop <- my_hits$Orthogroup
to_crop <- as.data.frame(to_crop)
colnames(to_crop) <- "Orthogroup"
my_hits <- my_hits[my_hits$picmin_fdr < 0.1,]
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


# NEW
my_res <- readRDS("orthogroup_results.rds")



for (og in my_hits$Orthogroup) {
og <- as.data.frame(og)
colnames(og) <- "Orthogroup"
 og_res <- og %>%
   left_join(my_res, by = c("Orthogroup"))
   og_res <-og_res[complete.cases(og_res), ]
   my_species <- unique(og_res$species)
   filtered_data <- my_res %>%
   group_by(Orthogroup) %>%
   filter(all(sort(my_species) %in% sort(species)) & length(species) == length(my_species) )
   matching_OGs <- unique(filtered_data$Orthogroup)
   matching_OGs <- as.data.frame(matching_OGs)
   colnames(matching_OGs) <- "Orthogroup"
   dist_dist <- matching_OGs %>%
     left_join(aa_dist, by = c("Orthogroup"))
	 
    hitter <- og %>%
    left_join(aa_dist, by = c("Orthogroup"))
	quantile_95 <- quantile(dist_dist$distance,0.95)
	
    filename <- paste0("plot_", og, ".svg")
	svg(filename, height =4, width =4)
     print(ggplot(dist_dist, aes(x = distance)) +
       geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +

       # Add vertical line at the mean
       geom_vline(data = hitter, aes(xintercept = distance), color = "red") +

       # Add vertical line at the 95th percentile
       geom_vline(xintercept = quantile_95, color = "black") +

       # Set labels and title
       labs(title = og,
            x = "Distance",
            y = "Frequency") + theme_bw())
  # Add your code here to perform operations on 'og'
  dev.off()
}
