# Read the Tau scores into my_data
# map is what I use to link TAIR_gene to Orthogroup_ID

my_data <- readRDS("Athal_tau_data.rds")
map <- readRDS("OG_map_Athal_COMBINED25species_updatedOF_221213.rds")


library(qvalue)

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



# Here below I simply join each row of my_data, which are A.thaliana genes, with the row with the same TAIR_gene in the map.
# I do this so that I am adding the Orthogroup ID to my_data, which I rename "all"
library(dplyr)

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
# the mean is not on 0 (but on -0.60) and if you take random draws of size 1000 (thats more or less the size of my singificant picmin results) from this Z score column,
# the mean keeps on staying around -0.60

final <- cbind(ortho_adjusted,z)
final <- final[is.finite(final$z),]





###DRAWS

my_list <- replicate(10000,sample_n(final, 1000))

final_to_plot <- c()


for (i in 1:10000) {
    final_to_plot <- append(final_to_plot,mean(my_list[,i]$z))
}
