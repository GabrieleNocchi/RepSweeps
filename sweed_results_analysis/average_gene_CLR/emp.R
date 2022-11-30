sweed <- read.table("final_genes_average_ortho.txt", header=FALSE)

sw <- sweed[,4]
ortho <- sweed[,5]


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


sweed <- subset(sweed, select = -c(V4,V5) )
sweed <- cbind(sweed,mean_emp_p, ortho)

write.table(sweed, file = "final_genes_average_ortho_tmp.txt", quote = FALSE, row.names = FALSE,col.names=FALSE, sep = "\t")

