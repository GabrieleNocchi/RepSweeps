library(dplyr)

genes <- read.table("final_genes_analysis.txt")

genes_minimum <- genes %>% group_by(V1) %>% slice(which.min(V2))

genes_minimum <- data.frame(genes_minimum)

genes_mean <- genes %>% group_by(V1) %>% summarise_at(vars(V2), list(name = mean))

genes_mean <- as.data.frame(genes_mean)

genes_adjusted <- cbind(genes_minimum, genes_mean$name)

colnames(genes_adjusted) <- c("gene","min_CLR", "scan_n", "mean_CLR")


# Og map
my_map  <- read.table("map.txt")
colnames(my_map) <- c("gene", "orthogroup")


# Retrieve OG for my genes
genes_with_OG <- genes_adjusted %>%
  left_join(my_map, by = c("gene"))

genes_with_OG <- genes_with_OG[complete.cases(genes_with_OG), ]


ortho_count <- table(genes_with_OG$orthogroup)

ortho_count <- as.data.frame(ortho_count)
colnames(ortho_count) <- c("orthogroup", "ortho_size")


genes_with_OG_paralogues_count <- genes_with_OG %>%
  left_join(ortho_count, by = c("orthogroup"))



# Remove large OGs
  genes_with_OG_paralogues_count <- genes_with_OG_paralogues_count[genes_with_OG_paralogues_count$ortho_size < 11,]
genes_mean_clr <- genes_with_OG_paralogues_count$mean_CLR

# Assign emp p values to average CLR per gene

  assign.pvalues <- function(array){
    #array <- sample(sw, 1000)
    pvalues <- array(0, length(array))

    ordered.indexes <- order(array)

    j <- length(array)
    for( i in ordered.indexes ){
      pvalues[i] <- j/length(array)
      j <- j-1
    }

    return(pvalues)
  }


  emp_p <- assign.pvalues(genes_mean_clr)


  genes_with_OG_paralogues_count$mean_CLR <- emp_p

  # mean_emp_p is the mean_CLR of the gene converted to emp_p value, hence mean_emp_p

  colnames(genes_with_OG_paralogues_count) <- c("gene", "min_CLR", "scan_n", "mean_emp_p", "orthogroup", "ortho_size" )
 


write.table(genes_with_OG_paralogues_count, file = "all_genes_average.txt", row.names = F, col.names = T, quote = F, sep = "\t")
