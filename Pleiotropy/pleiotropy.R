# Read the Tau scores into my_data
# map is what I use to link TAIR_gene to Orthogroup_ID
my_hits <- readRDS("gab_picmin_results.rds")
my_hits <- my_hits$picmin_res
my_hits <- my_hits[my_hits$picmin_fdr < 0.4,]
my_hits <- my_hits$Orthogroup
my_hits <- as.data.frame(my_hits)
colnames(my_hits) <- "Orthogroup"


map <- readRDS("OG_map_Athal_COMBINED25species_updatedOF_221213.rds")

assign.pvalues <- function(array){
  #array <- sample(sw, 1000)
  pvalues <- array(0, length(array))

  ordered.indexes <- order(array, decreasing = F)

  j <- length(array)
  for( i in ordered.indexes ){
    pvalues[i] <- j/length(array)
    j <- j-1
  }

  return(pvalues)
}

library(dplyr)
library(ggplot2)
library(qvalue)


############################# 1
my_data <- readRDS("Athal_tau_data.rds")


# Here below I simply join each row of my_data, which are A.thaliana genes, with the row with the same TAIR_gene in the map.
# I do this so that I am adding the Orthogroup ID to my_data, which I rename "all"


all <- my_data %>%
  left_join(map, by = c("TAIR_gene"))

# This line below as there are some NAs in the table that cause issues later
all<-all[complete.cases(all), ]

# Here I take tau_nolog column and transform it into emp-p
annot <- all[,4]

mean_tau_emp_p <- assign.pvalues(annot)

### last column of all now has empirical p values for tau_nolog

all <- cbind(all,mean_tau_emp_p)


# Here I simply count the size (number of paralogues) of each orthogroup, by simply checking the frequency of each unique orthogroup ID in "all"
ortho_count <- table(all$Orthogroup)

ortho_count <- as.data.frame(ortho_count)
colnames(ortho_count) <- c("Orthogroup", "Freq")

# each row of all, corresponding to genes, now gets another column also showing the orthogroup size count (Freq)

all <- all %>%
  left_join(ortho_count, by = c("Orthogroup"))


### DS Correction per Orthogroup
### I take only one gene/row per orthogroup, the row with the gene with minimum emp_p in that Orthogroup

ortho_minimum <- all %>% group_by(Orthogroup) %>% slice(which.min(mean_tau_emp_p))

ortho_minimum <- data.frame(ortho_minimum)

gabriele_dunnsidak <- function(x, y){
    1 - ((1-x)^y)
}


# I correct the minimum emp with DS based on the number of paralogues, which is stored in the column Freq

dunnsidak_orthopvalues <- gabriele_dunnsidak(ortho_minimum$mean_tau_emp_p,ortho_minimum$Freq)

ortho_adjusted <- cbind(ortho_minimum, dunnsidak_orthopvalues)

ortho_adjusted$dunnsidak_orthopvalues <- empPvals(-ortho_adjusted$dunnsidak_orthopvalues,-ortho_adjusted$dunnsidak_orthopvalues)

## Now I Transorm the DS adjusted ortho emp p to Z scores)
z <- qnorm(1 - ortho_adjusted$dunnsidak_orthopvalues)


# Final file. Last column has the Z scores for all the Orthogroup in A thaliana based on tau_nolog

final <- cbind(ortho_adjusted,z)
final <- final[is.finite(final$z),]

my_hits_z <- my_hits %>%
  left_join(final, by = c("Orthogroup"))
my_hits_z1 <-my_hits_z[complete.cases(my_hits_z), ]



###DRAWS

my_list <- replicate(1000,sample_n(final, 120))

final_to_plot <- c()


for (i in 1:1000) {
    final_to_plot <- append(final_to_plot,mean(my_list[,i]$z))
}


cat <- rep("Tissue Specificity - Tau", length(final_to_plot))

df <- cbind(final_to_plot,cat)

# Plot
df1 <- as.data.frame(df)




############################# 2
my_data <- readRDS("220525_Athal_coexpression_node_stats.rds")


# Here below I simply join each row of my_data, which are A.thaliana genes, with the row with the same TAIR_gene in the map.
# I do this so that I am adding the Orthogroup ID to my_data, which I rename "all"


all <- my_data %>%
  left_join(map, by = c("TAIR_gene"))

# This line below as there are some NAs in the table that cause issues later
all<-all[complete.cases(all), ]

# Here I take 2nd column and transform it into emp-p
annot <- all[,2]

mean_tau_emp_p <- assign.pvalues(annot)

### last column of all now has empirical p values for 2nd col

all <- cbind(all,mean_tau_emp_p)


# Here I simply count the size (number of paralogues) of each orthogroup, by simply checking the frequency of each unique orthogroup ID in "all"
ortho_count <- table(all$Orthogroup)

ortho_count <- as.data.frame(ortho_count)
colnames(ortho_count) <- c("Orthogroup", "Freq")

# each row of all, corresponding to genes, now gets another column also showing the orthogroup size count (Freq)

all <- all %>%
  left_join(ortho_count, by = c("Orthogroup"))


### DS Correction per Orthogroup
### I take only one gene/row per orthogroup, the row with the gene with minimum emp_p in that Orthogroup

ortho_minimum <- all %>% group_by(Orthogroup) %>% slice(which.min(mean_tau_emp_p))

ortho_minimum <- data.frame(ortho_minimum)

gabriele_dunnsidak <- function(x, y){
    1 - ((1-x)^y)
}


# I correct the minimum emp with DS based on the number of paralogues, which is stored in the column Freq

dunnsidak_orthopvalues <- gabriele_dunnsidak(ortho_minimum$mean_tau_emp_p,ortho_minimum$Freq)

ortho_adjusted <- cbind(ortho_minimum, dunnsidak_orthopvalues)

ortho_adjusted$dunnsidak_orthopvalues <- empPvals(-ortho_adjusted$dunnsidak_orthopvalues,-ortho_adjusted$dunnsidak_orthopvalues)

## Now I Transorm the DS adjusted ortho emp p to Z scores)
z <- qnorm(1 - ortho_adjusted$dunnsidak_orthopvalues)


# Final file. Last column has the Z scores for all the Orthogroup in A thaliana based on 2nd col

final <- cbind(ortho_adjusted,z)
final <- final[is.finite(final$z),]

my_hits_z <- my_hits %>%
  left_join(final, by = c("Orthogroup"))
my_hits_z2 <-my_hits_z[complete.cases(my_hits_z), ]



my_list <- replicate(1000,sample_n(final, 120))

final_to_plot <- c()


for (i in 1:1000) {
    final_to_plot <- append(final_to_plot,mean(my_list[,i]$z))
}


cat <- rep("Node Betweenness", length(final_to_plot))

df <- cbind(final_to_plot,cat)

# Plot
df2 <- as.data.frame(df)




############################# 3

my_data <- readRDS("220525_Athal_coexpression_node_stats.rds")


# Here below I simply join each row of my_data, which are A.thaliana genes, with the row with the same TAIR_gene in the map.
# I do this so that I am adding the Orthogroup ID to my_data, which I rename "all"


all <- my_data %>%
  left_join(map, by = c("TAIR_gene"))

# This line below as there are some NAs in the table that cause issues later
all<-all[complete.cases(all), ]

# Here I take 3rd col column and transform it into emp-p
annot <- all[,3]

mean_tau_emp_p <- assign.pvalues(annot)

### last column of all now has empirical p values for 3rd col

all <- cbind(all,mean_tau_emp_p)


# Here I simply count the size (number of paralogues) of each orthogroup, by simply checking the frequency of each unique orthogroup ID in "all"
ortho_count <- table(all$Orthogroup)

ortho_count <- as.data.frame(ortho_count)
colnames(ortho_count) <- c("Orthogroup", "Freq")

# each row of all, corresponding to genes, now gets another column also showing the orthogroup size count (Freq)

all <- all %>%
  left_join(ortho_count, by = c("Orthogroup"))


### DS Correction per Orthogroup
### I take only one gene/row per orthogroup, the row with the gene with minimum emp_p in that Orthogroup

ortho_minimum <- all %>% group_by(Orthogroup) %>% slice(which.min(mean_tau_emp_p))

ortho_minimum <- data.frame(ortho_minimum)

gabriele_dunnsidak <- function(x, y){
    1 - ((1-x)^y)
}


# I correct the minimum emp with DS based on the number of paralogues, which is stored in the column Freq

dunnsidak_orthopvalues <- gabriele_dunnsidak(ortho_minimum$mean_tau_emp_p,ortho_minimum$Freq)

ortho_adjusted <- cbind(ortho_minimum, dunnsidak_orthopvalues)

ortho_adjusted$dunnsidak_orthopvalues <- empPvals(-ortho_adjusted$dunnsidak_orthopvalues,-ortho_adjusted$dunnsidak_orthopvalues)

## Now I Transorm the DS adjusted ortho emp p to Z scores)
z <- qnorm(1 - ortho_adjusted$dunnsidak_orthopvalues)


# Final file. Last column has the Z scores for all the Orthogroup in A thaliana based on third col

final <- cbind(ortho_adjusted,z)
final <- final[is.finite(final$z),]

my_hits_z <- my_hits %>%
  left_join(final, by = c("Orthogroup"))
my_hits_z3 <-my_hits_z[complete.cases(my_hits_z), ]



my_list <- replicate(1000,sample_n(final, 120))

final_to_plot <- c()


for (i in 1:1000) {
    final_to_plot <- append(final_to_plot,mean(my_list[,i]$z))
}


cat <- rep("Node Strength", length(final_to_plot))

df <- cbind(final_to_plot,cat)

# Plot
df3 <- as.data.frame(df)




############################# 4
my_data <- readRDS("220525_Athal_coexpression_node_stats.rds")


# Here below I simply join each row of my_data, which are A.thaliana genes, with the row with the same TAIR_gene in the map.
# I do this so that I am adding the Orthogroup ID to my_data, which I rename "all"


all <- my_data %>%
left_join(map, by = c("TAIR_gene"))

# This line below as there are some NAs in the table that cause issues later
all<-all[complete.cases(all), ]

# Here I take 4th col column and transform it into emp-p
annot <- all[,4]

mean_tau_emp_p <- assign.pvalues(annot)

### last column of all now has empirical p values for 4th col

all <- cbind(all,mean_tau_emp_p)


# Here I simply count the size (number of paralogues) of each orthogroup, by simply checking the frequency of each unique orthogroup ID in "all"
ortho_count <- table(all$Orthogroup)

ortho_count <- as.data.frame(ortho_count)
colnames(ortho_count) <- c("Orthogroup", "Freq")

# each row of all, corresponding to genes, now gets another column also showing the orthogroup size count (Freq)

all <- all %>%
left_join(ortho_count, by = c("Orthogroup"))


### DS Correction per Orthogroup
### I take only one gene/row per orthogroup, the row with the gene with minimum emp_p in that Orthogroup

ortho_minimum <- all %>% group_by(Orthogroup) %>% slice(which.min(mean_tau_emp_p))

ortho_minimum <- data.frame(ortho_minimum)

gabriele_dunnsidak <- function(x, y){
1 - ((1-x)^y)
}


# I correct the minimum emp with DS based on the number of paralogues, which is stored in the column Freq

dunnsidak_orthopvalues <- gabriele_dunnsidak(ortho_minimum$mean_tau_emp_p,ortho_minimum$Freq)

ortho_adjusted <- cbind(ortho_minimum, dunnsidak_orthopvalues)

ortho_adjusted$dunnsidak_orthopvalues <- empPvals(-ortho_adjusted$dunnsidak_orthopvalues,-ortho_adjusted$dunnsidak_orthopvalues)

## Now I Transorm the DS adjusted ortho emp p to Z scores)
z <- qnorm(1 - ortho_adjusted$dunnsidak_orthopvalues)


# Final file. Last column has the Z scores for all the Orthogroup in A thaliana based on 4th col

final <- cbind(ortho_adjusted,z)
final <- final[is.finite(final$z),]

my_hits_z <- my_hits %>%
left_join(final, by = c("Orthogroup"))
my_hits_z4 <-my_hits_z[complete.cases(my_hits_z), ]



my_list <- replicate(1000,sample_n(final, 120))

final_to_plot <- c()


for (i in 1:1000) {
final_to_plot <- append(final_to_plot,mean(my_list[,i]$z))
}


cat <- rep("Node Degree", length(final_to_plot))

df <- cbind(final_to_plot,cat)

# Plot
df4 <- as.data.frame(df)




############################# 5
my_data <- readRDS("220525_Athal_coexpression_node_stats.rds")


# Here below I simply join each row of my_data, which are A.thaliana genes, with the row with the same TAIR_gene in the map.
# I do this so that I am adding the Orthogroup ID to my_data, which I rename "all"


all <- my_data %>%
  left_join(map, by = c("TAIR_gene"))

# This line below as there are some NAs in the table that cause issues later
all<-all[complete.cases(all), ]

# Here I take 5th col column and transform it into emp-p
annot <- all[,5]

mean_tau_emp_p <- assign.pvalues(annot)

### last column of all now has empirical p values for 5th col

all <- cbind(all,mean_tau_emp_p)


# Here I simply count the size (number of paralogues) of each orthogroup, by simply checking the frequency of each unique orthogroup ID in "all"
ortho_count <- table(all$Orthogroup)

ortho_count <- as.data.frame(ortho_count)
colnames(ortho_count) <- c("Orthogroup", "Freq")

# each row of all, corresponding to genes, now gets another column also showing the orthogroup size count (Freq)

all <- all %>%
  left_join(ortho_count, by = c("Orthogroup"))


### DS Correction per Orthogroup
### I take only one gene/row per orthogroup, the row with the gene with minimum emp_p in that Orthogroup

ortho_minimum <- all %>% group_by(Orthogroup) %>% slice(which.min(mean_tau_emp_p))

ortho_minimum <- data.frame(ortho_minimum)

gabriele_dunnsidak <- function(x, y){
    1 - ((1-x)^y)
}


# I correct the minimum emp with DS based on the number of paralogues, which is stored in the column Freq

dunnsidak_orthopvalues <- gabriele_dunnsidak(ortho_minimum$mean_tau_emp_p,ortho_minimum$Freq)

ortho_adjusted <- cbind(ortho_minimum, dunnsidak_orthopvalues)

ortho_adjusted$dunnsidak_orthopvalues <- empPvals(-ortho_adjusted$dunnsidak_orthopvalues,-ortho_adjusted$dunnsidak_orthopvalues)
## Now I Transorm the DS adjusted ortho emp p to Z scores)
z <- qnorm(1 - ortho_adjusted$dunnsidak_orthopvalues)


# Final file. Last column has the Z scores for all the Orthogroup in A thaliana based on 5th col

final <- cbind(ortho_adjusted,z)
final <- final[is.finite(final$z),]

my_hits_z <- my_hits %>%
  left_join(final, by = c("Orthogroup"))
my_hits_z5 <-my_hits_z[complete.cases(my_hits_z), ]



my_list <- replicate(1000,sample_n(final, 120))

final_to_plot <- c()


for (i in 1:1000) {
    final_to_plot <- append(final_to_plot,mean(my_list[,i]$z))
}


cat <- rep("Node Closeness", length(final_to_plot))

df <- cbind(final_to_plot,cat)







# Plot All
df5 <- as.data.frame(df)
df <- rbind(df1,df2,df3,df4,df5)

df$final_to_plot <- as.numeric(df$final_to_plot)
df1$final_to_plot <- as.numeric(df1$final_to_plot)
df2$final_to_plot <- as.numeric(df2$final_to_plot)
df3$final_to_plot <- as.numeric(df3$final_to_plot)
df4$final_to_plot <- as.numeric(df4$final_to_plot)
df5$final_to_plot <- as.numeric(df5$final_to_plot)



    ggplot(df, aes(cat, mean(final_to_plot))) +        # ggplot2 plot with confidence intervals
    geom_point(fill="black", color="black", size=4, shape = 16) +
    geom_linerange(aes(x = "Tissue Specificity - Tau",ymin = quantile(df1$final_to_plot,0.05), ymax = quantile(df1$final_to_plot,0.95)), size = 1.5) +
    geom_linerange(aes(x = "Node Betweenness",ymin = quantile(df2$final_to_plot,0.05), ymax = quantile(df2$final_to_plot,0.95)), size = 1.5) +
    geom_linerange(aes(x = "Node Strength",ymin = quantile(df3$final_to_plot,0.05), ymax = quantile(df3$final_to_plot,0.95)), size = 1.5) +
    geom_linerange(aes(x = "Node Degree",ymin = quantile(df4$final_to_plot,0.05), ymax = quantile(df4$final_to_plot,0.95)), size = 1.5) +
    geom_linerange(aes(x = "Node Closeness",ymin = quantile(df5$final_to_plot,0.05), ymax = quantile(df5$final_to_plot,0.95)), size = 1.5) +
    theme_classic() + coord_flip() +
    geom_point(aes(y=mean(my_hits_z1$z), x = "Tissue Specificity - Tau"),colour="red", shape = 17, size = 4) +
    geom_point(aes(y=mean(my_hits_z2$z), x = "Node Betweenness"),colour="red", shape = 17, size = 4) +
    geom_point(aes(y=mean(my_hits_z3$z), x = "Node Strength"),colour="red", shape = 17, size = 4) +
    geom_point(aes(y=mean(my_hits_z4$z), x = "Node Degree"),colour="red", shape = 17, size = 4) +
    geom_point(aes(y=mean(my_hits_z5$z), x = "Node Closeness"),colour="red", shape = 17, size = 4) +
    theme(axis.title.y=element_blank(),axis.line.y=element_blank(),
          axis.ticks.y=element_blank()) + ylab("Mean Orthogroup Z score") + ylim(-0.3,0.3) + geom_hline(yintercept=0, linetype="dashed") + geom_hline(yintercept=c(-0.2,-0.15,-0.1,0.1,0.15,0.2), linetype="dotted")
