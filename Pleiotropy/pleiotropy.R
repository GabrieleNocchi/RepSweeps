my_data <- readRDS("Athal_tau_data.rds")
map <- readRDS("OG_map_Athal_220927.rds")


annot <- my_data[,2]


assign.pvalues <- function(array){
  #array <- sample(sw, 1000)
  pvalues <- array(0, length(array))

  ordered.indexes <- order(array, decreasing = TRUE)

  j <- length(array)
  for( i in ordered.indexes ){
    pvalues[i] <- j/length(array)
    j <- j-1
  }

  return(pvalues)
}


mean_tau_emp_p <- assign.pvalues(annot)



my_data <- cbind(my_data,mean_tau_emp_p)


library(dplyr)

all <- my_data %>%
  left_join(map, by = c("TAIR_gene"))


all<-all[complete.cases(all), ]

ortho_count <- table(all$Orthogroup)

ortho_count <- as.data.frame(ortho_count)
colnames(ortho_count) <- c("Orthogroup", "Freq")

all <- all %>%
  left_join(ortho_count, by = c("Orthogroup"))




### DS Correction per Orthogroup


ortho_minimum <- all %>% group_by(Orthogroup) %>% slice(which.min(mean_tau_emp_p))

ortho_minimum <- data.frame(ortho_minimum)

gabriele_dunnsidak <- function(x, y){
    1 - ((1-x)^y)
}


dunnsidak_orthopvalues <- gabriele_dunnsidak(ortho_minimum$mean_tau_emp_p,ortho_minimum$Freq)

ortho_adjusted <- cbind(ortho_minimum, dunnsidak_orthopvalues)

z <- qnorm(ortho_adjusted$dunnsidak_orthopvalues)


final <- cbind(ortho_adjusted,z)
final <- final[is.finite(final$z),]
