sweed <- read.table("tmp_file.txt", header=FALSE)

sweed <- sweed[sweed$V6 < 11,]
sw <- sweed[,4]



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


sweed[,4] <- mean_emp_p 
write.table(sweed, file = "tmp_file.txt", quote = FALSE, row.names = FALSE,col.names=FALSE, sep = "\t")

