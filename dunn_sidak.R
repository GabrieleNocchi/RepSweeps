library(dplyr)

genes <- read.table("final_genes_analysis.txt")

genes_minimum <- genes %>% group_by(V1) %>% slice(which.min(V2))

genes_minimum <- data.frame(genes_minimum)

gabriele_dunnsidak <- function(x, y){
    1 - ((1-x)^y)
}


adjusted_pvalues <- gabriele_dunnsidak(genes_minimum$V2,genes_minimum$V3)


genes_adjusted <- cbind(genes_minimum, adjusted_pvalues)
colnames(genes_adjusted) <- c("gene","min_emp_p", "scan_n", "dunn_sidak")


write.table(genes_adjusted, file = "final_genes_dunak.txt", quote = FALSE, row.names = FALSE,col.names=TRUE, sep = "\t")
