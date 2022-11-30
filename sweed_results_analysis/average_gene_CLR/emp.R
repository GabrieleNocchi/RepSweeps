sweed <- read.table("final_genes_average.txt", header=TRUE)

sw <- sweed$mean_CLR



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


mean_emp_p <- assign.pvalues(sw)

sweed <- cbind(sweed,mean_emp_p)
sweed <- subset(sweed, select = -c(mean_CLR) )

write.table(sweed, file = "final_genes_average_2.txt", quote = FALSE, row.names = FALSE,col.names=FALSE, sep = "\t")

