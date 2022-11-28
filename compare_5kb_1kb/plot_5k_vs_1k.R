setwd("/users/gnocc/Desktop")
library(dplyr)

df <- read.table("to_plot_in_R.txt")
df2 <- df %>% group_by(V1) %>% summarise_at(vars(V6), list(name = mean))

df3 <- df %>% group_by(V1) %>% slice(which.max(V6))


a <- as.data.frame(df2)
b <- as.data.frame(df3)
c <- cbind(b,a)

colnames(c) <- c("window_ID1", "start", "end", "null","CLR_window", "maximum_CLR_window", "window_ID2", "mean_CLR")

plot(-log10(b$V5), -log10(b$V6), xlab = "-log10 5kb CLR", ylab = "-log10 1kb_max_CLR")
cor(b$V5,b$V6)

plot(-log10(c$CLR_window), -log10(c$mean_CLR), xlab = "-log10 5kb CLR", ylab = "-log10 1kb_average_CLR")
cor(c$CLR_window,c$mean_CLR)
