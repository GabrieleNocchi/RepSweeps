library(RColorBrewer)
library(ggplot2)

a <- read.table("unrelated_to_plot.txt", sep = "\t", h = T)


mycolors = c("lightgoldenrod1", "darkgrey")



ggplot(a, aes(fill=cat, y=count, x=genome)) +
  geom_bar(position="stack", stat="identity") + coord_flip() + theme_classic() + xlab("Genome") + ylab("Individuals Count") +
    scale_fill_viridis_d(option="cividis") + labs(fill='Individuals') + ggtitle("Number of selected individuals per dataset kin < 0.2")
