library(ggplot2)
data <-readRDS("gab_picmin_results.rds")
data <- data$picmin_res


fir <- data[data$picmin_fdr <= 0.5 & data$picmin_fdr > 0.4,]
sec <- data[data$picmin_fdr <= 0.4 & data$picmin_fdr > 0.3,]
thi <- data[data$picmin_fdr <= 0.3 & data$picmin_fdr > 0.2,]
fou <- data[data$picmin_fdr <= 0.2 & data$picmin_fdr > 0.1,]
fif <- data[data$picmin_fdr <= 0.1,]

a <- length(fir$p)
b <- length(sec$p)
c <- length(thi$p)
d <- length(fou$p)
e <- length(fif$p)


count <- c(a,b,c,d,e)
FDR <- c("\u2264 0.5","\u2264 0.4", "\u2264 0.3", "\u2264 0.2", "\u2264 0.1")
picmin <- rep("PicMin", 5)


my_data <- data.frame(count, FDR, picmin)

p1 <- ggplot(my_data, aes(fill=FDR, y=count, x=picmin)) +
      geom_bar(position="stack", stat="identity", width = 0.5) + theme_classic() +
      theme(axis.title.y=element_blank(),
        axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + ylab("OG count\n") +
        scale_colour_viridis_d(direction = -1) + scale_fill_viridis_d(direction = -1, option  = "cividis") +
        coord_flip() +
        scale_y_continuous(breaks=seq(0,500,20), position = "right") +
        ggtitle(paste("PicMin OGs = ", format(round(as.numeric(length(data$p)), 1), big.mark=",")))



library(ggplot2)
library(RColorBrewer)
library(RColorBrewer)
library(colorspace)

pal <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF", "#393B79", "#637939", "#8C6D31", "gold", "#4A7EBB", "#AEC7E8", "#FFBB78")


pic <- readRDS("gab_picmin_results.rds")
pic <- pic$picmin_res
pic <- pic[pic$picmin_fdr < 0.4,]
results <- readRDS("orthogroup_results.rds")


library(dplyr)
data <- pic %>% left_join(results, by = c("Orthogroup"))


p2<- ggplot(data, aes(x=reorder(species, ortho_DS, FUN = median), y=ortho_DS, fill = species)) +
     geom_boxplot() + theme_classic()  + theme(legend.position = "none") + ylab("Ortho_DS\n") + ggtitle("Distribution of the OGs DS corrected emp p-values across species for the PicMin OGs with FDR < 0.4 (280 OGs)") + xlab("Species") +
     scale_fill_manual(values=pal)



library(dplyr)
low_p_count <- data %>% group_by(species) %>% count(ortho_DS < 0.1)
low_p_count <- data.frame(low_p_count)
colnames(low_p_count) <- c("species", "count", "n")

below <- low_p_count[low_p_count$count == "TRUE",]
above <- low_p_count[low_p_count$count == "FALSE",]

below <- cbind (below$species, below$n)
colnames(below) <- c("species", "below")

above <- cbind(above$species, above$n)
colnames(above) <- c("species", "above")

above <- as.data.frame(above)
below <- as.data.frame(below)

reformatted_data <- below %>% left_join(above, by = c("species"))
contribution <- as.numeric(reformatted_data$below) / (as.numeric(reformatted_data$below)+ as.numeric(reformatted_data$above))
final_data <- cbind(reformatted_data,contribution)


# p3 <- ggplot(final_data, aes(x = reorder(species, contribution),y = 1, fill = contribution)) + scale_fill_viridis_c(direction = 1, option  = "cividis") + theme_classic() +
#   geom_tile(height =0.5) +
#   theme(axis.title.y=element_blank(),
#         axis.line.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank()) +
#         xlab("Species") + labs(fill='Species contribution') + ggtitle("Species contribution to PicMin OGs with FDR < 0.4 --> N(OG-p < 0.1)/N(OG-p)")


p4 <- ggplot(data=final_data, aes(x=reorder(species,contribution), y=contribution*100)) +
geom_bar(stat="identity", width = 0.6,fill = "gold", col = "black") + coord_flip() + theme_classic() + xlab("Species") + ylab("Contribution %") + labs(fill='Species') + theme(legend.position = "none") +
ggtitle("Species contribution to PicMin OGs with FDR < 0.4 --> N(OG-p < 0.1)/N(OG-p)")


library(gridExtra)
grid.arrange(p1, p4,p2,heights = c(1,1, 1))
