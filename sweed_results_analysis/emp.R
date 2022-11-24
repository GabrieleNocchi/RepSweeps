sweed <- read.table("all.bed", header=FALSE)

sw <- sweed[,5]



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


sw.pval <- assign.pvalues(sw)

sweed <- cbind(sweed,sw.pval)


write.table(sweed, file = "all_emp.bed", quote = FALSE, row.names = FALSE,col.names=FALSE, sep = "\t")
