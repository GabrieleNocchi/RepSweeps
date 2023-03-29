assign.pvalues <- function(array){
  #array <- sample(sw, 1000)
  pvalues <- array(0, length(array))

  ordered.indexes <- order(array, decreasing = TRUE)

  j <- length(array)
  for( i in ordered.indexes ){
    pvalues[i] <- j/length(array)
    j <- j-1
  }

  return(pvalues)
}


a <- read.table("kubota_ahalleri_final_orthogroup.txt", h = T)
species <- rep("ahalleri", length(a$orthogroup))
a <- cbind(a,species)
pvals <- a$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
a$dunnsidak_orthopvalues <- tmp_emp
write.table(a, file = "kubota_ahalleri_final_orthogroup.txt", quote = FALSE, row.names = FALSE,col.names=TRUE, sep = "\t")


b <- read.table("ingvarsson_ptremula_final_orthogroup.txt", h = T)
species <- rep("ptremula", length(b$orthogroup))
b <- cbind(b,species)
pvals <- b$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
b$dunnsidak_orthopvalues <- tmp_emp
write.table(b, file = "ingvarsson_ptremula_final_orthogroup.txt", quote = FALSE, row.names = FALSE,col.names=TRUE, sep = "\t")


c <- read.table("lowry_phallii_final_orthogroup.txt", h = T)
species <- rep("phalli", length(c$orthogroup))
c <- cbind(c,species)
pvals <- c$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
c$dunnsidak_orthopvalues <- tmp_emp
write.table(c, file = "lowry_phallii_final_orthogroup.txt", quote = FALSE, row.names = FALSE,col.names=TRUE, sep = "\t")


d <- read.table("mitchell_bstricta_final_orthogroup.txt", h = T)
species <- rep("bstricta", length(d$orthogroup))
d <- cbind(d,species)
pvals <- d$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
d$dunnsidak_orthopvalues <- tmp_emp
write.table(d, file = "mitchell_bstricta_final_orthogroup.txt", quote = FALSE, row.names = FALSE,col.names=TRUE, sep = "\t")


e <- read.table("murray_ealb_final_orthogroup.txt", h = T)
species <- rep("ealb", length(e$orthogroup))
e <- cbind(e,species)
pvals <- e$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
e$dunnsidak_orthopvalues <- tmp_emp
write.table(e, file = "murray_ealb_final_orthogroup.txt", quote = FALSE, row.names = FALSE,col.names=TRUE, sep = "\t")


f <- read.table("murray_esid_final_orthogroup.txt", h = T)
species <- rep("esid", length(f$orthogroup))
f <- cbind(f,species)
pvals <- f$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
f$dunnsidak_orthopvalues <- tmp_emp
write.table(f, file = "murray_esid_final_orthogroup.txt", quote = FALSE, row.names = FALSE,col.names=TRUE, sep = "\t")


g <- read.table("murray_emag_final_orthogroup.txt", h = T)
species <- rep("emag", length(g$orthogroup))
g <- cbind(g,species)
pvals <- g$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
g$dunnsidak_orthopvalues <- tmp_emp
write.table(g, file = "murray_emag_final_orthogroup.txt", quote = FALSE, row.names = FALSE,col.names=TRUE, sep = "\t")


h <- read.table("todesco_hpet_final_orthogroup.txt", h = T)
species <- rep("hpet", length(h$orthogroup))
h <- cbind(h,species)
pvals <- h$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
h$dunnsidak_orthopvalues <- tmp_emp
write.table(h, file = "todesco_hpet_final_orthogroup.txt", quote = FALSE, row.names = FALSE,col.names=TRUE, sep = "\t")


i <- read.table("tiffin_mtruncatula_final_orthogroup.txt", h = T)
species <- rep("mtruncatula", length(i$orthogroup))
i <- cbind(i,species)
pvals <- i$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
i$dunnsidak_orthopvalues <- tmp_emp
write.table(i, file = "tiffin_mtruncatula_final_orthogroup.txt", quote = FALSE, row.names = FALSE,col.names=TRUE, sep = "\t")


j <- read.table("weigel_athaliana_IBE_final_orthogroup.txt", h = T)
species <- rep("athaliana", length(j$orthogroup))
j <- cbind(j,species)
pvals <- j$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
j$dunnsidak_orthopvalues <- tmp_emp
write.table(j, file = "weigel_athaliana_IBE_final_orthogroup.txt", quote = FALSE, row.names = FALSE,col.names=TRUE, sep = "\t")


k <- read.table("weigel_capsella_final_orthogroup.txt", h = T)
species <- rep("capsella", length(k$orthogroup))
k <- cbind(k,species)
pvals <- k$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
k$dunnsidak_orthopvalues <- tmp_emp
write.table(k, file = "weigel_capsella_final_orthogroup.txt", quote = FALSE, row.names = FALSE,col.names=TRUE, sep = "\t")


l <- read.table("wright_atuber_final_orthogroup.txt", h = T)
species <- rep("atuber", length(l$orthogroup))
l <- cbind(l,species)
pvals <- l$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
l$dunnsidak_orthopvalues <- tmp_emp
write.table(l, file = "wright_atuber_final_orthogroup.txt", quote = FALSE, row.names = FALSE,col.names=TRUE, sep = "\t")


m <- read.table("todesco_hann_final_orthogroup.txt", h = T)
species <- rep("hann", length(m$orthogroup))
m <- cbind(m,species)
pvals <- m$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
m$dunnsidak_orthopvalues <- tmp_emp
write.table(m, file = "todesco_hann_final_orthogroup.txt", quote = FALSE, row.names = FALSE,col.names=TRUE, sep = "\t")


n <- read.table("todesco_hargo_final_orthogroup.txt", h = T)
species <- rep("hargo", length(n$orthogroup))
n <- cbind(n,species)
pvals <- n$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
n$dunnsidak_orthopvalues <- tmp_emp
write.table(n, file = "todesco_hargo_final_orthogroup.txt", quote = FALSE, row.names = FALSE,col.names=TRUE, sep = "\t")


o <- read.table("evans_ptricho_final_orthogroup.txt", h = T)
species <- rep("ptricho", length(o$orthogroup))
o <- cbind(o,species)
pvals <- o$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
o$dunnsidak_orthopvalues <- tmp_emp
write.table(o, file = "evans_ptricho_final_orthogroup.txt", quote = FALSE, row.names = FALSE,col.names=TRUE, sep = "\t")


p <- read.table("nocchi_bplaty_final_orthogroup.txt", h = T)
species <- rep("bplaty", length(p$orthogroup))
p <- cbind(p,species)
pvals <- p$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
p$dunnsidak_orthopvalues <- tmp_emp
write.table(p, file = "nocchi_bplaty_final_orthogroup.txt", quote = FALSE, row.names = FALSE,col.names=TRUE, sep = "\t")


q <- read.table("salojarvi_bpendula_final_orthogroup.txt", h = T)
species <- rep("bpendula", length(q$orthogroup))
q <- cbind(q,species)
pvals <- q$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
q$dunnsidak_orthopvalues <- tmp_emp
write.table(q, file = "salojarvi_bpendula_final_orthogroup.txt", quote = FALSE, row.names = FALSE,col.names=TRUE, sep = "\t")


results <- rbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q)
climate <- rep("mean_temp", length(results$species))

results <- cbind(results,climate)
colnames(results) <- c("gene","min_CLR","scan_n","mean_emp_p","Orthogroup","ortho_size", "ortho_DS", "species", "climate")
saveRDS(results, "orthogroup_results.rds")
