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



b <- read.table("ingvarsson_ptremula_final_orthogroup.txt", h = T)
species <- rep("ptremula", length(b$orthogroup))
b <- cbind(b,species)
pvals <- b$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
b$dunnsidak_orthopvalues <- tmp_emp



c <- read.table("lowry_phallii_final_orthogroup.txt", h = T)
species <- rep("phalli", length(c$orthogroup))
c <- cbind(c,species)
pvals <- c$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
c$dunnsidak_orthopvalues <- tmp_emp



d <- read.table("mitchell_bstricta_final_orthogroup.txt", h = T)
species <- rep("bstricta", length(d$orthogroup))
d <- cbind(d,species)
pvals <- d$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
d$dunnsidak_orthopvalues <- tmp_emp



e <- read.table("murray_ealb_final_orthogroup.txt", h = T)
species <- rep("ealb", length(e$orthogroup))
e <- cbind(e,species)
pvals <- e$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
e$dunnsidak_orthopvalues <- tmp_emp



f <- read.table("murray_esid_final_orthogroup.txt", h = T)
species <- rep("esid", length(f$orthogroup))
f <- cbind(f,species)
pvals <- f$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
f$dunnsidak_orthopvalues <- tmp_emp



g <- read.table("murray_emag_final_orthogroup.txt", h = T)
species <- rep("emag", length(g$orthogroup))
g <- cbind(g,species)
pvals <- g$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
g$dunnsidak_orthopvalues <- tmp_emp



h <- read.table("todesco_hpet_final_orthogroup.txt", h = T)
species <- rep("hpet", length(h$orthogroup))
h <- cbind(h,species)
pvals <- h$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
h$dunnsidak_orthopvalues <- tmp_emp



i <- read.table("tiffin_mtruncatula_final_orthogroup.txt", h = T)
species <- rep("mtruncatula", length(i$orthogroup))
i <- cbind(i,species)
pvals <- i$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
i$dunnsidak_orthopvalues <- tmp_emp



j <- read.table("weigel_athaliana_IBE_final_orthogroup.txt", h = T)
species <- rep("athaliana", length(j$orthogroup))
j <- cbind(j,species)
pvals <- j$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
j$dunnsidak_orthopvalues <- tmp_emp



k <- read.table("weigel_capsella_final_orthogroup.txt", h = T)
species <- rep("capsella", length(k$orthogroup))
k <- cbind(k,species)
pvals <- k$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
k$dunnsidak_orthopvalues <- tmp_emp



l <- read.table("wright_atuber_final_orthogroup.txt", h = T)
species <- rep("atuber", length(l$orthogroup))
l <- cbind(l,species)
pvals <- l$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
l$dunnsidak_orthopvalues <- tmp_emp



m <- read.table("todesco_hann_final_orthogroup.txt", h = T)
species <- rep("hann", length(m$orthogroup))
m <- cbind(m,species)
pvals <- m$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
m$dunnsidak_orthopvalues <- tmp_emp



n <- read.table("todesco_hargo_final_orthogroup.txt", h = T)
species <- rep("hargo", length(n$orthogroup))
n <- cbind(n,species)
pvals <- n$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
n$dunnsidak_orthopvalues <- tmp_emp


results <- rbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n)
climate <- rep("mean_temp", length(results$species))

results <- cbind(results,climate)
colnames(results) <- c("gene","min_CLR","scan_n","mean_emp_p","Orthogroup","ortho_size", "ortho_DS", "species", "climate")
saveRDS(results, "orthogroup_results.rds")
