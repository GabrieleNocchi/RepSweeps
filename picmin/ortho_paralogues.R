library(dplyr)
library(ggplot2)

results <- readRDS("orthogroup_results.rds")
results$species <- sub("ptricho", "P. trichocarpa", results$species)
results$species <- sub("emag", "E. magnificata", results$species)
results$species <- sub("ealb", "E. albens", results$species)
results$species <- sub("esid", "E. sideroxylon", results$species)
results$species <- sub("bplaty", "B. platyphylla", results$species)
results$species <- sub("bpendula", "B. pendula", results$species)
results$species <- sub("hann", "H. annus", results$species)
results$species <- sub("hargo", "H. argophyllus", results$species)
results$species <- sub("hpet", "H. petiolaris", results$species)
results$species <- sub("capsella", "C. rubella", results$species)
results$species <- sub("ahalleri", "A. halleri", results$species)
results$species <- sub("athaliana", "A. thaliana", results$species)
results$species <- sub("ptremula", "P. tremula", results$species)
results$species <- sub("bstricta", "B. stricta", results$species)
results$species <- sub("phalli", "P. hallii", results$species)
results$species <- sub("atuber", "A. tuberculatus", results$species)
results$species <- sub("mtruncatula", "M. truncatula", results$species)

my_data <- results %>% group_by(species, ortho_size) %>% summarise(total_count=n())

my_data <- as.data.frame(my_data)



p1 <-   ggplot(my_data, aes(fill=as.factor(ortho_size), y=total_count, x=species)) +
        geom_bar(stat="identity", width = 0.8, position = position_stack(reverse = TRUE)) + coord_flip() + theme_classic() + scale_fill_viridis_d(direction = -1, option  = "mako") +
        labs(fill='Paralogues N') + ylab("Number") + xlab("Species") + ggtitle("Paralogues per orthogroup") +
        theme(legend.key.size = unit(0.2, 'cm'), plot.title = element_text(size = 10, face = "bold"))




my_data_2 <- results %>% group_by(Orthogroup) %>% summarise(total_count=n())


my_data_2 <- data.frame(my_data_2)

total <- table(my_data_2$total_count)

total <- as.data.frame(total)

p2 <- ggplot(data=total, aes(x=Var1, y=Freq, fill = "constant")) +
      geom_bar(stat="identity", color= "black", width = 0.8,linewidth = 0.3) + theme_classic() + scale_fill_viridis_d(direction = -1, option  = "mako") + xlab("Number of Species") + ylab("Number") + ggtitle("OGs tested in PicMin")+
      theme(plot.title = element_text(size = 10, face = "bold"), legend.position = "none")

ortho_species <- results %>% left_join(my_data_2, by = c("Orthogroup"))

my_data_3 <- ortho_species %>% group_by(species, total_count) %>% summarise(number=n())

my_data_3 <- as.data.frame(my_data_3)


p3 <- ggplot(my_data_3, aes(fill=as.factor(total_count), y=number, x=species)) +
      geom_bar(stat="identity", width = 0.8) + coord_flip() + theme_classic() + scale_fill_viridis_d(direction = 1, option  = "mako") +
      labs(fill='Species N') + ylab("Number") + xlab("Species") + ggtitle("Species per orthogroup") +
      theme(legend.key.size = unit(0.2, 'cm'),plot.title = element_text(size = 10, face = "bold"))




svg("OGs_assement.svg", height = 7, width = 7)
library(gridExtra)
grid.arrange(p1,p3,p2,ncol=2)
dev.off()
