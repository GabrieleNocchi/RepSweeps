a <- readRDS("gab_picmin_results.rds")
a <- a$picmin_res
#a <- a[a$picmin_fdr < 0.5,]
a <- a[,1]

my_results <- readRDS("orthogroup_results.rds")

for (i in 1:length(a)) {
    selection <- my_results[my_results$Orthogroup == a[i],]
    genes <- selection$gene
    genes <- unique(sort(genes))
    write.table(genes,file = paste0(a[i], ".txt"), quote = FALSE, row.names = F, col.names = F)
}
