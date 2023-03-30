# Read the Tau scores into my_data
# map is what I use to link TAIR_gene to Orthogroup_ID
my_hits <- readRDS("gab_picmin_results.rds")
my_hits <- my_hits$picmin_res
# my_hits <- my_hits[my_hits$picmin_fdr < 0.5,]
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


z1 <- sum(my_hits_z1$z)/(sqrt(length(my_hits_z1$z)))




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




z2 <- sum(my_hits_z2$z)/(sqrt(length(my_hits_z2$z)))






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




z3 <- sum(my_hits_z3$z)/(sqrt(length(my_hits_z3$z)))





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



z4 <- sum(my_hits_z4$z)/(sqrt(length(my_hits_z4$z)))





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




z5 <- sum(my_hits_z5$z)/(sqrt(length(my_hits_z5$z)))







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



z6 <- sum(my_hits_z6$z)/(sqrt(length(my_hits_z6$z)))






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




z7 <- sum(my_hits_z7$z)/(sqrt(length(my_hits_z7$z)))





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



z8 <- sum(my_hits_z8$z)/(sqrt(length(my_hits_z8$z)))




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

z9 <- sum(my_hits_z9$z)/(sqrt(length(my_hits_z9$z)))
#
# par(mfrow=c(3,3))
# plot(-log(my_hits_z1$picmin_fdr),my_hits_z1$z, xlab = "-log PicMin FDR", ylab = "Tau Z")
# abline(h=0, col = "red")
# plot(-log(my_hits_z2$picmin_fdr),my_hits_z2$z, xlab = "-log PicMin FDR", ylab = "Athaliana betweenness Z")
# abline(h=0, col = "red")
# plot(-log(my_hits_z3$picmin_fdr),my_hits_z3$z, xlab = "-log PicMin FDR", ylab = "Athaliana strength Z")
# abline(h=0, col = "red")
# plot(-log(my_hits_z4$picmin_fdr),my_hits_z4$z, xlab = "-log PicMin FDR", ylab = "Athaliana degree Z")
# abline(h=0, col = "red")
# plot(-log(my_hits_z5$picmin_fdr),my_hits_z5$z, xlab = "-log PicMin FDR", ylab = "Athaliana closeness Z")
# abline(h=0, col = "red")
# plot(-log(my_hits_z6$picmin_fdr),my_hits_z6$z, xlab = "-log PicMin FDR", ylab = "Medicago betweenness Z")
# abline(h=0, col = "red")
# plot(-log(my_hits_z7$picmin_fdr),my_hits_z7$z, xlab = "-log PicMin FDR", ylab = "Medicago stength Z")
# abline(h=0, col = "red")
# plot(-log(my_hits_z8$picmin_fdr),my_hits_z8$z, xlab = "-log PicMin FDR", ylab = "Medicago degree Z")
# abline(h=0, col = "red")
# plot(-log(my_hits_z9$picmin_fdr),my_hits_z9$z, xlab = "-log PicMin FDR", ylab = "Medicago closeness Z")
# abline(h=0, col = "red")



library(ggplot2)
library(gridExtra)

# Define function to create individual plots
create_plot <- function(data, xlab, ylab) {
  ggplot(data, aes(x = -log(picmin_fdr), y = z)) +
    geom_point( color = "black", fill= "lightgoldenrod", shape = 21, size = 2) +
    geom_hline(yintercept = 0, color = "red") +
    xlab(xlab) +
    ylab(ylab) +
    theme_classic()
}

# Create individual plots
plot1 <- create_plot(my_hits_z1, "-log PicMin FDR", "Tau Z")
plot2 <- create_plot(my_hits_z2, "-log PicMin FDR", "Athaliana betweenness Z")
plot3 <- create_plot(my_hits_z3, "-log PicMin FDR", "Athaliana strength Z")
plot4 <- create_plot(my_hits_z4, "-log PicMin FDR", "Athaliana degree Z")
plot5 <- create_plot(my_hits_z5, "-log PicMin FDR", "Athaliana closeness Z")
plot6 <- create_plot(my_hits_z6, "-log PicMin FDR", "Medicago betweenness Z")
plot7 <- create_plot(my_hits_z7, "-log PicMin FDR", "Medicago strength Z")
plot8 <- create_plot(my_hits_z8, "-log PicMin FDR", "Medicago degree Z")
plot9 <- create_plot(my_hits_z9, "-log PicMin FDR", "Medicago closeness Z")

# Arrange plots in 3x3 grid using gridExtra
grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, ncol = 3)
