library(ggplot2)
library(reshape2)
library(dplyr)
library(circlize)


# HITS
z <- readRDS("gab_picmin_results.rds")
z <- z$picmin_res
z <- z[z$picmin_fdr<0.2,]
z <- as.data.frame(cbind(z$Orthogroup, z$climate_var))
colnames(z) <- c("Orthogroup", "removeme")


df <- readRDS("orthogroup_results.rds")
df <- df %>%   left_join(z, by = c("Orthogroup"))
df <- df[complete.cases(df),]

df$difference <- sapply(strsplit(sub(".*:", "", df$gene), "-"), function(x) as.numeric(x[2]) - as.numeric(x[1]))


hits_mean <- df %>% group_by(Orthogroup) %>% summarise_at(vars(difference), list(name = mean))
hits_mean <- as.data.frame(hits_mean)


# ALL
z <- readRDS("gab_picmin_results.rds")
z <- z$picmin_res
z <- as.data.frame(cbind(z$Orthogroup, z$climate_var))
colnames(z) <- c("Orthogroup", "removeme")


df <- readRDS("orthogroup_results.rds")
df <- df %>%   left_join(z, by = c("Orthogroup"))
df <- df[complete.cases(df),]



df$difference <- sapply(strsplit(sub(".*:", "", df$gene), "-"), function(x) as.numeric(x[2]) - as.numeric(x[1]))

genes_mean <- df %>% group_by(Orthogroup) %>% summarise_at(vars(difference), list(name = mean))
genes_mean <- as.data.frame(genes_mean)


### DRAWS

my_list <- replicate(10000,sample_n(genes_mean, length(hits_mean$Orthogroup)))

final_to_plot <- c()


for (i in 1:10000) {
    final_to_plot <- append(final_to_plot,mean(my_list[,i]$name))
}





cat <- rep("Genes Length", length(final_to_plot))

df <- cbind(final_to_plot,cat)

# Plot
df <- as.data.frame(df)

df$final_to_plot <- as.numeric(df$final_to_plot)





a <- mean(as.numeric(df$final_to_plot))
a1 <- quantile(df$final_to_plot,0.025)
a2 <- quantile(df$final_to_plot,0.975)


df <- cbind(a,"Genes Length")
colnames(df) <- c("final_to_plot", "cat")


df <-as.data.frame(df)
df$final_to_plot <- as.numeric(df$final_to_plot)



library(ggplot2)
ggplot(df, aes(cat, (final_to_plot))) +        # ggplot2 plot with confidence intervals
    geom_point(color = "black", fill= "lightgoldenrod", shape = 19, size = 4) +
    geom_linerange(aes(x = "Genes Length",ymin = a1, ymax = a2), size = 1.5) +
	 theme_classic() + coord_flip() +
    geom_point(aes(y=mean(hits_mean$name), x = "Genes Length"),colour="black", fill = "red", shape = 24, size = 3) +
    theme(axis.title.y=element_blank(),axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.title = element_text(hjust = 0.5, vjust = 1.5)) +
    ylab("Mean genes length in OGs") +  ggtitle("FDR < 0.2")
