library(data.table)
library(dplyr)
library(ggplot2)
# Step 1: Read the results file
results <- fread("fst_results.weir.fst", header=TRUE)

# Step 2: Read the GFF file and extract gene coordinates
gff <- read.table("genes.gff", h =T, col.names=c("CHROM", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))

# Extract gene coordinates
genes <- gff[, c("CHROM", "start", "end")]

# Step 3: Calculate average FST for each gene
gene_avg_FST <- numeric(nrow(genes))

for (i in seq_along(genes$CHROM)) {
  chrom <- genes$CHROM[i]
  start <- genes$start[i]
  end <- genes$end[i]

  # Subset results for the current chromosome and within gene coordinates
  subset_results <- results[results$CHROM == chrom & results$POS >= start & results$POS <= end, ]

  # Calculate average FST
  avg_FST <- mean(subset_results$WEIR_AND_COCKERHAM_FST, na.rm=TRUE)
  gene_avg_FST[i] <- avg_FST
  print(avg_FST)
  print(i)
}

# Step 4: Combine gene information with average FST
final_output <- cbind(genes, mean_FST=gene_avg_FST)

# Print or save final_output as needed
final_output$gene <- sprintf("%s:%s-%s", final_output$CHROM, final_output$start, final_output$end)
final_output <- final_output %>%
  mutate(mean_FST = ifelse(mean_FST < 0, 0, mean_FST))
final_output <- final_output[complete.cases(final_output),]

OG <-  readRDS("gab_picmin_results.rds")
OG <- OG$picmin_res
OG <- OG[OG$picmin_fdr < 0.1,]
my_hits <- data.frame()

my_results <- readRDS("orthogroup_results.rds")
for (row_id in OG$Orthogroup) {

selection <- my_results[my_results$Orthogroup == row_id & my_results$ortho_DS < 0.1 & my_results$species == "bpendula",]
my_hits <- rbind(my_hits, selection)
}


print(my_hits)
head(final_output)


my_hits <- my_hits %>%
  left_join(final_output, by = c("gene"))
my_hits<-my_hits[complete.cases(my_hits), ]

pdf("platyphylla_fst.pdf", height = 3, width = 5)
ggplot(final_output, aes(x = mean_FST)) +
geom_histogram(fill = "lightgoldenrod", color = "black", binwidth=0.005) + geom_vline(xintercept = my_hits$mean_FST, color= "darkred", size = 0.5) + theme_bw() + geom_vline(xintercept = quantile(final_output$mean_FST, 0.95), color= "black", size = 0.5) +
xlab("FST across genes") + ylab("Frequency")
dev.off()



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


emp_p <- assign.pvalues(final_output$mean_FST)

final_output <- cbind(final_output, emp_p)

my_hits <- my_hits %>%
  left_join(final_output, by = c("gene"))
my_hits<-my_hits[complete.cases(my_hits), ]
write.table(file = "platyphylla_fst.txt", my_hits, row.names = F, quote = F)

