library(dplyr)

df <- read.table("to_plot_in_R.txt")
df2 <- df %>% group_by(V1) %>% summarise_at(vars(V7), list(name = mean))

df3 <- df %>% group_by(V1) %>% slice(which.max(V7))


a <- as.data.frame(df2)
b <- as.data.frame(df3)
c <- cbind(b,a)

colnames(c) <- c("window_ID1", "start", "end", "null","CLR_window", "SNP_count", "maximum_CLR_window", "window_ID2", "mean_CLR")


plot(log10(c$CLR_window), log10(c$mean_CLR), xlab = "log10 5kb CLR", ylab = "log10 1kb_average_CLR")
cor(c$CLR_window,c$mean_CLR)



library(dplyr)
library(ggplot2)

# c <- c %>%
#   filter(CLR_window > 0 & mean_CLR > 0)


p <- ggplot(c, aes(log10(CLR_window + 0.00001), log10(mean_CLR + 0.00001)))
p + geom_point(shape =21, size = 3, color = "black", fill = "lightgoldenrod") + theme_classic() + xlab("log10 10kb windows (5kb spaced scans) CLR") + ylab("log10 1kb-calculated average CLR per 10kb window")
