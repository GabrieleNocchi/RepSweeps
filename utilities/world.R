setwd("/Users/gnocc/Desktop/RepSweeps/world_map/data/")
library("rnaturalearth")
library("rnaturalearthdata")
library(RColorBrewer)
library(data.table)
library(ggplot2)
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)





    list_csv_files <- list.files(path = "/Users/gnocc/Desktop/RepSweeps/world_map/data/")
    df2 = do.call(rbind, lapply(list_csv_files, function(x) fread(x, stringsAsFactors = FALSE, h = F, sep = ",",select=c(1,2,3,4,5,6,7))))
    df2

mycolors = c("#1B9E77","#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D" ,"#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "cadetblue", "#FB9A99", "#E31A1C")

ggplot(data = world) +
    geom_sf() + theme_classic() + geom_point(data = df2, aes(x = V7, y = V6, col = V1), size = 0.5,
        shape = 19) + scale_color_manual(values = mycolors)+ guides(colour = guide_legend(override.aes = list(size=5))) +
        theme(axis.title.y=element_blank(),
        axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.line.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), plot.title = element_text(hjust = 0.5)) + labs(col='Species') + ggtitle("RepSweeps WGS Datasets")
