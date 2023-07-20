library(ape)
tree<-read.tree("alignment.nwk")
my_m <-cophenetic(tree)

my_m <-as.data.frame.table(my_m)

write.table(my_m, file = "formatted_phylogenies.txt", row.names = FALSE, quote = FALSE)
