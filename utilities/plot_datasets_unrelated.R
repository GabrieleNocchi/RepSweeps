library(RColorBrewer)
library(ggplot2)

a <- read.table("unrelated_to_plot.txt", sep = "\t", h = T)


mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6))



ggplot(a, aes(fill=cat, y=count, x=genome)) +
  geom_bar(position="stack", stat="identity") + coord_flip() + theme_classic() + xlab("Genome") + ylab("Individuals Count") +
    scale_fill_manual(values = mycolors) + labs(fill='Individuals') + ggtitle("Number of selected individuals per dataset kin < 0.2") 
