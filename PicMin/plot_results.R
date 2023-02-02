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

ggplot(my_data, aes(fill=FDR, y=count, x=picmin)) +
    geom_bar(position="stack", stat="identity", width = 0.25) + theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + ylab("OG count\n") +
        scale_colour_viridis_d(direction = -1) + scale_fill_viridis_d(direction = -1) +
        coord_flip() +
        scale_y_continuous(breaks=seq(0,280,20), position = "right") +
        ggtitle(paste("Total number of OGs tested in PicMin = ", format(round(as.numeric(length(data$p)), 1), big.mark=",")))
