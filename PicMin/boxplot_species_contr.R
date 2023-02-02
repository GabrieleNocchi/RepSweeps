library(ggplot2)
library(RColorBrewer)

pic <- readRDS("gab_picmin_results.rds")
pic <- pic$picmin_res
pic <- pic[pic$picmin_fdr < 0.5,]
results <- readRDS("orthogroup_results.rds")


library(dplyr)
data <- pic %>% left_join(results, by = c("Orthogroup"))
mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6))

ggplot(data, aes(x=reorder(species, ortho_DS, FUN = median), y=ortho_DS, fill = species)) +
    geom_boxplot(width = 0.5) + theme_classic()  + theme(legend.position = "none") + ylab("Ortho_DS\n") + ggtitle("Distribution of the OGs DS corrected emp p-values across species for the PicMin OGs with FDR < 0.5 (280 OGs)") + xlab("Species") +
    scale_fill_manual(values = mycolors)
