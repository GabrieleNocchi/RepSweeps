library(dplyr)
library(ggplot2)

results <- readRDS("orthogroup_results.rds")

my_data <- results %>% group_by(species, ortho_size) %>% summarise(total_count=n())

my_data <- as.data.frame(my_data)



p1 <-   ggplot(my_data, aes(fill=as.factor(ortho_size), y=total_count, x=species)) +
        geom_bar(stat="identity", width = 0.8, position = position_stack(reverse = TRUE)) + coord_flip() + theme_classic() + scale_fill_viridis_d(direction = -1, option  = "mako") +
        labs(fill='Paralogues N') + ylab("Number") + xlab("Species") + ggtitle("Paralogues Number per Orthogroup") +
        theme(legend.key.size = unit(0.2, 'cm'), plot.title = element_text(size = 10, face = "bold"))




my_data_2 <- results %>% group_by(Orthogroup) %>% summarise(total_count=n())


my_data_2 <- data.frame(my_data_2)

total <- table(my_data_2$total_count)

total <- as.data.frame(total)

p2 <- ggplot(data=total, aes(x=Var1, y=Freq, fill = Var1)) +
      geom_bar(stat="identity") + coord_flip() + theme_classic() + scale_fill_viridis_d(direction = -1, option  = "mako") + xlab("Number of Orthologues") + ylab("Number") + ggtitle("Total Number of OGs by orthologues Number\nPicMin used only OGs with at least 7 Orthologues")+
      labs(fill='Orthologues N') +
      theme(legend.key.size = unit(0.2, 'cm'),plot.title = element_text(size = 10, face = "bold"))

ortho_species <- results %>% left_join(my_data_2, by = c("Orthogroup"))

my_data_3 <- ortho_species %>% group_by(species, total_count) %>% summarise(number=n())

my_data_3 <- as.data.frame(my_data_3)


p3 <- ggplot(my_data_3, aes(fill=as.factor(total_count), y=number, x=species)) +
      geom_bar(stat="identity", width = 0.8) + coord_flip() + theme_classic() + scale_fill_viridis_d(direction = 1, option  = "mako") +
      labs(fill='Species N') + ylab("Number") + xlab("Species") + ggtitle("Species Number per Orthogroup in each Species") +
      theme(legend.key.size = unit(0.2, 'cm'),plot.title = element_text(size = 10, face = "bold"))




svg("OGs_assement.svg")
library(gridExtra)
grid.arrange(p1,p3,p2,ncol=1,heights = c(2, 2, 2))
dev.off()
