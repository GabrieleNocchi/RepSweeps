# Read the Tau scores into my_data
# map is what I use to link TAIR_gene to Orthogroup_ID
my_hits <- readRDS("gab_picmin_results.rds")
my_hits <- my_hits$picmin_res
my_hits <- my_hits[my_hits$picmin_fdr < 0.1,]
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







my_data_to_plot <- rbind(z1,z2,z3,z4,z5,z6,z7,z8,z9)
my_categories <- c("Arabidopsis Tissue Specificity - Tau","Arabidopsis Node Betweenness","Arabidopsis Node Strength","Arabidopsis Node Degree","Arabidopsis Node Closeness","Medicago Node Betweenness","Medicago Node Strength","Medicago Node Degree","Medicago Node Closeness")

my_data_to_plot <- cbind(my_data_to_plot,my_categories)
colnames(my_data_to_plot) <- c("stouf","cat")
my_data_to_plot <- as.data.frame(my_data_to_plot)
my_data_to_plot$stouf <- as.numeric(my_data_to_plot$stouf)

ggplot(my_data_to_plot, aes(y=stouf, x=cat)) +
    geom_bar(stat="identity", color = "black", fill = "lightgoldenrod") + theme_classic() + coord_flip() + geom_hline(yintercept=1.95, linetype="dashed") + geom_hline(yintercept=-1.95, linetype="dashed") +
    scale_y_continuous(breaks = c(-4,-3,-2,-1,1,2,3,4), limit = c(-5,5)) + ylab("Stouffer's Z") + theme(axis.title.y=element_blank(),axis.line.y=element_blank(),axis.ticks.y=element_blank(), plot.title = element_text(hjust = 0.5, vjust = 3.5)) + theme(aspect.ratio = .4) +
    annotate("text",x="" ,y=3.5,label="High Pleiotropy",fontface = "bold") + annotate("text",x="",y=-3.5,label="Low Pleiotropy",fontface = "bold") +  ggtitle("FDR < 0.1")
