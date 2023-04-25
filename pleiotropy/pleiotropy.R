# Read the Tau scores into my_data
# map is what I use to link TAIR_gene to Orthogroup_ID
my_hits <- readRDS("gab_picmin_results.rds")
my_hits <- my_hits$picmin_res
my_hits <- my_hits[my_hits$picmin_fdr < 0.5,]
my_hits <- my_hits$Orthogroup
my_hits <- as.data.frame(my_hits)
colnames(my_hits) <- "Orthogroup"


to_crop <- readRDS("gab_picmin_results.rds")
to_crop <- to_crop$picmin_res

map <- read.table("Athal_map_OG_GAB.txt", fill = TRUE, h = T)
colnames(map) <- c("ID","GENE_COORD","TAIR_gene","GENE_NAME", "Orthogroup")

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


assign.pvalues_2 <- function(array){
  #array <- sample(sw, 1000)
  pvalues <- array(0, length(array))

  ordered.indexes <- order(array, decreasing = T)

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

#### Taking care of using only the OGs I tested in PicMin
all <- all %>%
  left_join(to_crop, by = c("Orthogroup"))
  all<-all[complete.cases(all), ]
  all <- all[,1:15]

# Here I take tau_nolog column and transform it into emp-p
annot <- all[,4]

mean_tau_emp_p <- assign.pvalues_2(annot)

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



### DRAWS

my_list <- replicate(10000,sample_n(final, length(my_hits_z1$Orthogroup)))

final_to_plot <- c()


for (i in 1:10000) {
    final_to_plot <- append(final_to_plot,mean(my_list[,i]$z))
}


cat <- rep("Arabidopsis Tissue Specificity - Tau", length(final_to_plot))

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


#### Taking care of using only the OGs I tested in PicMin
all <- all %>%
  left_join(to_crop, by = c("Orthogroup"))
  all<-all[complete.cases(all), ]
  all <- all[,1:14]

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



### DRAWS

my_list <- replicate(10000,sample_n(final, length(my_hits_z2$Orthogroup)))

final_to_plot <- c()


for (i in 1:10000) {
    final_to_plot <- append(final_to_plot,mean(my_list[,i]$z))
}


cat <- rep("Arabidopsis Node Betweenness", length(final_to_plot))

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

#### Taking care of using only the OGs I tested in PicMin
all <- all %>%
  left_join(to_crop, by = c("Orthogroup"))
  all<-all[complete.cases(all), ]
  all <- all[,1:14]


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



### DRAWS

my_list <- replicate(10000,sample_n(final, length(my_hits_z3$Orthogroup)))

final_to_plot <- c()


for (i in 1:10000) {
    final_to_plot <- append(final_to_plot,mean(my_list[,i]$z))
}


cat <- rep("Arabidopsis Node Strength", length(final_to_plot))

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


#### Taking care of using only the OGs I tested in PicMin
all <- all %>%
  left_join(to_crop, by = c("Orthogroup"))
  all<-all[complete.cases(all), ]
  all <- all[,1:14]



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



### DRAWS

my_list <- replicate(10000,sample_n(final, length(my_hits_z4$Orthogroup)))

final_to_plot <- c()


for (i in 1:10000) {
final_to_plot <- append(final_to_plot,mean(my_list[,i]$z))
}


cat <- rep("Arabidopsis Node Degree", length(final_to_plot))

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

#### Taking care of using only the OGs I tested in PicMin
all <- all %>%
  left_join(to_crop, by = c("Orthogroup"))
  all<-all[complete.cases(all), ]
  all <- all[,1:14]

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



### DRAWS

my_list <- replicate(10000,sample_n(final, length(my_hits_z5$Orthogroup)))

final_to_plot <- c()


for (i in 1:10000) {
    final_to_plot <- append(final_to_plot,mean(my_list[,i]$z))
}


cat <- rep("Arabidopsis Node Closeness", length(final_to_plot))

df <- cbind(final_to_plot,cat)
df5 <- as.data.frame(df)






### MEDICAGO
map <- read.table("Mtrunc_map_OG_GAB.txt", fill = TRUE, h = T)
colnames(map) <- c("ID","GENE_COORD","GENE_NAME", "biomart_gene", "Orthogroup")
############################# 6
my_data <- readRDS("220927_Mtrunc_coexpression_node_stats.rds")


# Here below I simply join each row of my_data, which are A.thaliana genes, with the row with the same TAIR_gene in the map.
# I do this so that I am adding the Orthogroup ID to my_data, which I rename "all"


all <- my_data %>%
  left_join(map, by = c("biomart_gene"))

# This line below as there are some NAs in the table that cause issues later
all<-all[complete.cases(all), ]

#### Taking care of using only the OGs I tested in PicMin
all <- all %>%
  left_join(to_crop, by = c("Orthogroup"))
  all<-all[complete.cases(all), ]
  all <- all[,1:14]

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
my_hits_z6 <-my_hits_z[complete.cases(my_hits_z), ]



### DRAWS

my_list <- replicate(10000,sample_n(final, length(my_hits_z6$Orthogroup)))

final_to_plot <- c()


for (i in 1:10000) {
    final_to_plot <- append(final_to_plot,mean(my_list[,i]$z))
}


cat <- rep("Medicago Node Betweenness", length(final_to_plot))

df <- cbind(final_to_plot,cat)

# Plot
df6 <- as.data.frame(df)




############################# 7

my_data <- readRDS("220927_Mtrunc_coexpression_node_stats.rds")


# Here below I simply join each row of my_data, which are A.thaliana genes, with the row with the same TAIR_gene in the map.
# I do this so that I am adding the Orthogroup ID to my_data, which I rename "all"


all <- my_data %>%
  left_join(map, by = c("biomart_gene"))

# This line below as there are some NAs in the table that cause issues later
all<-all[complete.cases(all), ]

#### Taking care of using only the OGs I tested in PicMin
all <- all %>%
  left_join(to_crop, by = c("Orthogroup"))
  all<-all[complete.cases(all), ]
  all <- all[,1:14]

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
my_hits_z7 <-my_hits_z[complete.cases(my_hits_z), ]



### DRAWS

my_list <- replicate(10000,sample_n(final, length(my_hits_z7$Orthogroup)))

final_to_plot <- c()


for (i in 1:10000) {
    final_to_plot <- append(final_to_plot,mean(my_list[,i]$z))
}


cat <- rep("Medicago Node Strength", length(final_to_plot))

df <- cbind(final_to_plot,cat)

# Plot
df7 <- as.data.frame(df)




############################# 8
my_data <- readRDS("220927_Mtrunc_coexpression_node_stats.rds")


# Here below I simply join each row of my_data, which are A.thaliana genes, with the row with the same TAIR_gene in the map.
# I do this so that I am adding the Orthogroup ID to my_data, which I rename "all"


all <- my_data %>%
left_join(map, by = c("biomart_gene"))

# This line below as there are some NAs in the table that cause issues later
all<-all[complete.cases(all), ]

#### Taking care of using only the OGs I tested in PicMin
all <- all %>%
  left_join(to_crop, by = c("Orthogroup"))
  all<-all[complete.cases(all), ]
  all <- all[,1:14]

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
my_hits_z8 <-my_hits_z[complete.cases(my_hits_z), ]



### DRAWS

my_list <- replicate(10000,sample_n(final, length(my_hits_z8$Orthogroup)))

final_to_plot <- c()


for (i in 1:10000) {
final_to_plot <- append(final_to_plot,mean(my_list[,i]$z))
}


cat <- rep("Medicago Node Degree", length(final_to_plot))

df <- cbind(final_to_plot,cat)

# Plot
df8 <- as.data.frame(df)




############################# 9
my_data <- readRDS("220927_Mtrunc_coexpression_node_stats.rds")


# Here below I simply join each row of my_data, which are A.thaliana genes, with the row with the same TAIR_gene in the map.
# I do this so that I am adding the Orthogroup ID to my_data, which I rename "all"


all <- my_data %>%
  left_join(map, by = c("biomart_gene"))

# This line below as there are some NAs in the table that cause issues later
all<-all[complete.cases(all), ]

#### Taking care of using only the OGs I tested in PicMin
all <- all %>%
  left_join(to_crop, by = c("Orthogroup"))
  all<-all[complete.cases(all), ]
  all <- all[,1:14]

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
my_hits_z9 <-my_hits_z[complete.cases(my_hits_z), ]



### DRAWS

my_list <- replicate(10000,sample_n(final, length(my_hits_z9$Orthogroup)))

final_to_plot <- c()


for (i in 1:10000) {
    final_to_plot <- append(final_to_plot,mean(my_list[,i]$z))
}


cat <- rep("Medicago Node Closeness", length(final_to_plot))

df <- cbind(final_to_plot,cat)
df9 <- as.data.frame(df)

# Plot All

df <- rbind(df1,df2,df3,df4,df5,df6,df7,df8,df9)

df$final_to_plot <- as.numeric(df$final_to_plot)
df1$final_to_plot <- as.numeric(df1$final_to_plot)
df2$final_to_plot <- as.numeric(df2$final_to_plot)
df3$final_to_plot <- as.numeric(df3$final_to_plot)
df4$final_to_plot <- as.numeric(df4$final_to_plot)
df5$final_to_plot <- as.numeric(df5$final_to_plot)
df6$final_to_plot <- as.numeric(df6$final_to_plot)
df7$final_to_plot <- as.numeric(df7$final_to_plot)
df8$final_to_plot <- as.numeric(df8$final_to_plot)
df9$final_to_plot <- as.numeric(df9$final_to_plot)



### test to make plot less heavy for PDF
a <- mean(as.numeric(df1$final_to_plot))
a1 <- quantile(df1$final_to_plot,0.025)
a2 <- quantile(df1$final_to_plot,0.975)
b <- mean(as.numeric(df2$final_to_plot))
b1 <- quantile(df2$final_to_plot,0.025)
b2 <- quantile(df2$final_to_plot,0.975)
c <- mean(as.numeric(df3$final_to_plot))
c1 <- quantile(df3$final_to_plot,0.025)
c2 <- quantile(df3$final_to_plot,0.975)
d <- mean(as.numeric(df4$final_to_plot))
d1 <- quantile(df4$final_to_plot,0.025)
d2 <- quantile(df4$final_to_plot,0.975)
e <- mean(as.numeric(df5$final_to_plot))
e1 <- quantile(df5$final_to_plot,0.025)
e2 <- quantile(df5$final_to_plot,0.975)
f <- mean(as.numeric(df6$final_to_plot))
f1 <- quantile(df6$final_to_plot,0.025)
f2 <- quantile(df6$final_to_plot,0.975)
g <- mean(as.numeric(df7$final_to_plot))
g1 <- quantile(df7$final_to_plot,0.025)
g2 <- quantile(df7$final_to_plot,0.975)
h <- mean(as.numeric(df8$final_to_plot))
h1 <- quantile(df8$final_to_plot,0.025)
h2 <- quantile(df8$final_to_plot,0.975)
i <- mean(as.numeric(df9$final_to_plot))
i1 <- quantile(df9$final_to_plot,0.025)
i2 <- quantile(df9$final_to_plot,0.975)


df <- rbind(a,b,c,d,e,f,g,h,i)
cat <- c("Arabidopsis Tissue Specificity - Tau","Arabidopsis Node Betweenness","Arabidopsis Node Strength","Arabidopsis Node Degree","Arabidopsis Node Closeness","Medicago Node Betweenness","Medicago Node Strength","Medicago Node Degree","Medicago Node Closeness")
df <- cbind(df,cat)
colnames(df) <- c("final_to_plot", "cat")


df <-as.data.frame(df)
df$final_to_plot <- as.numeric(df$final_to_plot)


    ggplot(df, aes(cat, (final_to_plot))) +        # ggplot2 plot with confidence intervals
    geom_point(color = "black", fill= "lightgoldenrod", shape = 19, size = 4) +
    geom_linerange(aes(x = "Arabidopsis Tissue Specificity - Tau",ymin = a1, ymax = a2), size = 1.5) +
    geom_linerange(aes(x = "Arabidopsis Node Betweenness",ymin = b1, ymax = b2), size = 1.5) +
    geom_linerange(aes(x = "Arabidopsis Node Strength",ymin = c1, ymax = c2), size = 1.5) +
    geom_linerange(aes(x = "Arabidopsis Node Degree",ymin = d1, ymax = d2), size = 1.5) +
    geom_linerange(aes(x = "Arabidopsis Node Closeness",ymin = e1, ymax = e2), size = 1.5) +

    geom_linerange(aes(x = "Medicago Node Betweenness",ymin = f1, ymax = f2), size = 1.5) +
    geom_linerange(aes(x = "Medicago Node Strength",ymin = g1, ymax = g2), size = 1.5) +
    geom_linerange(aes(x = "Medicago Node Degree",ymin = h1, ymax = h2), size = 1.5) +
    geom_linerange(aes(x = "Medicago Node Closeness",ymin = i1, ymax = i2), size = 1.5) +
    theme_classic() + coord_flip() +
    geom_point(aes(y=mean(my_hits_z1$z), x = "Arabidopsis Tissue Specificity - Tau"),colour="black", fill = "red", shape = 24, size = 3) +
    geom_point(aes(y=mean(my_hits_z2$z), x = "Arabidopsis Node Betweenness"),colour="black", fill = "red",shape = 24, size = 3) +
    geom_point(aes(y=mean(my_hits_z3$z), x = "Arabidopsis Node Strength"),colour="black", fill = "red", shape = 24, size = 3) +
    geom_point(aes(y=mean(my_hits_z4$z), x = "Arabidopsis Node Degree"),colour="black", fill = "red",shape = 24, size = 3) +
    geom_point(aes(y=mean(my_hits_z5$z), x = "Arabidopsis Node Closeness"),colour="black", fill = "red",shape = 24, size = 3) +

    geom_point(aes(y=mean(my_hits_z6$z), x = "Medicago Node Betweenness"),colour="black", fill = "red",shape = 24, size = 3) +
    geom_point(aes(y=mean(my_hits_z7$z), x = "Medicago Node Strength"),colour="black", fill = "red",shape = 24, size = 3) +
    geom_point(aes(y=mean(my_hits_z8$z), x = "Medicago Node Degree"),colour="black", fill = "red", shape = 24, size = 3) +
    geom_point(aes(y=mean(my_hits_z9$z), x = "Medicago Node Closeness"),colour="black", fill = "red", shape = 24, size = 3) +
    theme(axis.title.y=element_blank(),axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.title = element_text(hjust = 0.5, vjust = 1.5)) +
    ylab("Mean Orthogroup Z score") + geom_hline(yintercept=0, linetype="dashed") + geom_hline(yintercept=c(-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2,0.25,0.3), linetype="dotted") +
    scale_y_continuous(breaks = c(-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6), limit = c(-0.6,0.6)) + annotate("text",x="" ,y=0.55,label="High Pleiotropy",fontface = "bold") + annotate("text",x="",y=-0.55,label="Low Pleiotropy",fontface = "bold") +  ggtitle("FDR < 0.5")
