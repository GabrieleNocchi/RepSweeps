library(dplyr)

# Taking average CLR for each gene, based on all the scans within that gene

genes <- read.table("final_genes_all_rec.txt")



genes_mean <- genes %>% group_by(V1) %>% summarise_at(vars(V2), list(name = mean))

genes_mean <- as.data.frame(genes_mean)


colnames(genes_mean) <- c("gene", "rec_mean")


# Og map
my_map  <- read.table("map.txt")
colnames(my_map) <- c("gene", "orthogroup")


# Retrieve OG for my genes
genes_with_OG <- genes_mean %>%
  left_join(my_map, by = c("gene"))

genes_with_OG <- genes_with_OG[complete.cases(genes_with_OG), ]



ortho_count <- table(genes_with_OG$orthogroup)

ortho_count <- as.data.frame(ortho_count)
colnames(ortho_count) <- c("orthogroup", "freq")


genes_with_OG_paralogues_count <- genes_with_OG %>%
  left_join(ortho_count, by = c("orthogroup"))



# Remove large OGs
  genes_with_OG_paralogues_count <- genes_with_OG_paralogues_count[genes_with_OG_paralogues_count$freq < 11,]
  genes_mean_rec <- genes_with_OG_paralogues_count$rec_mean




# Assign emp p values to average CLR per gene

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


  emp_p <- assign.pvalues(genes_mean_rec)


  genes_with_OG_paralogues_count$rec_mean <- emp_p

write.table(file = "final_rec_emp.txt",genes_with_OG_paralogues_count, col.names = TRUE, quote = FALSE, row.names = FALSE, sep = "\t")

my_results <- read.table("all_genes_average.txt", h = T)


comb <- genes_with_OG_paralogues_count %>% left_join(my_results, by = c("gene"))


library(ggdensity)

pdf("hargo.pdf")
ggplot(comb,aes(rec_mean,mean_emp_p)) +
    ggdensity::geom_hdr(
    aes(fill = after_stat(probs)), 
    alpha = 1) +
    theme_minimal()+
    theme(panel.grid = element_blank(),
          strip.text.y = element_text(angle = 0),
          strip.text.x = element_text(angle = 30)) +
    labs(y = "OmegaPlus ep-value",x = "Recombination ep-value")
dev.off()
