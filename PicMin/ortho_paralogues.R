library(dplyr)
library(ggplot2)
results <- readRDS("orthogroup_results.rds")

my_data <- results %>% group_by(species, ortho_size) %>%
    summarise(total_count=n())

    my_data <- as.data.frame(my_data)




    ggplot(my_data, aes(fill=as.factor(ortho_size), y=total_count, x=species)) +
    geom_bar(stat="identity", width = 0.8, position = position_stack(reverse = TRUE)) + coord_flip() + theme_classic() + scale_fill_viridis_d(direction = -1) +
    labs(fill='Paralogues N') + ylab("Number") + xlab("Species") + ggtitle("Paralogues Number per Orthogroup")
