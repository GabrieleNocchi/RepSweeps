library(dplyr)

ortho <- read.table("final_genes_dunak_ortho.txt", header = FALSE)


ortho_gene_count <- ortho %>% count(V5)

write.table(ortho_gene_count, file = "ortho_gene_count.txt", quote = FALSE, row.names = FALSE,col.names=FALSE, sep = "\t")
