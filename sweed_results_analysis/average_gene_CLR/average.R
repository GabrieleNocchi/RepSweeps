library(dplyr)

# Taking average CLR for each gene, based on all the scans within that gene

genes <- read.table("final_genes_analysis.txt")

genes_minimum <- genes %>% group_by(V1) %>% slice(which.min(V2))

genes_minimum <- data.frame(genes_minimum)

genes_mean <- genes %>% group_by(V1) %>% summarise_at(vars(V2), list(name = mean))

genes_mean <- as.data.frame(genes_mean)

genes_adjusted <- cbind(genes_minimum, genes_mean$name)

colnames(genes_adjusted) <- c("gene","min_emp_p", "scan_n", "mean_emp_p")


write.table(genes_adjusted, file = "final_genes_average.txt", quote = FALSE, row.names = FALSE,col.names=TRUE, sep = "\t")
