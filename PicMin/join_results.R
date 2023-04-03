tested <- read.table("tested_OGs.txt", h = T)
library(dplyr)
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
a <- tested %>%
  left_join(a, by = c("orthogroup"))
a<-a[complete.cases(a), ]
pvals <- a$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
a$dunnsidak_orthopvalues <- tmp_emp


b <- read.table("ingvarsson_ptremula_final_orthogroup.txt", h = T)
species <- rep("ptremula", length(b$orthogroup))
b <- cbind(b,species)
b <- tested %>%
  left_join(b, by = c("orthogroup"))
b<-b[complete.cases(b), ]
pvals <- b$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
b$dunnsidak_orthopvalues <- tmp_emp



c <- read.table("lowry_phallii_final_orthogroup.txt", h = T)
species <- rep("phalli", length(c$orthogroup))
c <- cbind(c,species)
c <- tested %>%
  left_join(c, by = c("orthogroup"))
c<-c[complete.cases(c), ]
pvals <- c$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
c$dunnsidak_orthopvalues <- tmp_emp



d <- read.table("mitchell_bstricta_final_orthogroup.txt", h = T)
species <- rep("bstricta", length(d$orthogroup))
d <- cbind(d,species)
d <- tested %>%
  left_join(d, by = c("orthogroup"))
d<-d[complete.cases(d), ]
pvals <- d$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
d$dunnsidak_orthopvalues <- tmp_emp



e <- read.table("murray_ealb_final_orthogroup.txt", h = T)
species <- rep("ealb", length(e$orthogroup))
e <- cbind(e,species)
e <- tested %>%
  left_join(e, by = c("orthogroup"))
e<-e[complete.cases(e), ]
pvals <- e$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
e$dunnsidak_orthopvalues <- tmp_emp



f <- read.table("murray_esid_final_orthogroup.txt", h = T)
species <- rep("esid", length(f$orthogroup))
f <- cbind(f,species)
f <- tested %>%
  left_join(f, by = c("orthogroup"))
f<-f[complete.cases(f), ]
pvals <- f$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
f$dunnsidak_orthopvalues <- tmp_emp



g <- read.table("murray_emag_final_orthogroup.txt", h = T)
species <- rep("emag", length(g$orthogroup))
g <- cbind(g,species)
g <- tested %>%
  left_join(g, by = c("orthogroup"))
g<-g[complete.cases(g), ]
pvals <- g$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
g$dunnsidak_orthopvalues <- tmp_emp



h <- read.table("todesco_hpet_final_orthogroup.txt", h = T)
species <- rep("hpet", length(h$orthogroup))
h <- cbind(h,species)
h <- tested %>%
  left_join(h, by = c("orthogroup"))
h<-h[complete.cases(h), ]
pvals <- h$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
h$dunnsidak_orthopvalues <- tmp_emp


i <- read.table("tiffin_mtruncatula_final_orthogroup.txt", h = T)
species <- rep("mtruncatula", length(i$orthogroup))
i <- cbind(i,species)
i <- tested %>%
  left_join(i, by = c("orthogroup"))
i<-i[complete.cases(i), ]
pvals <- i$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
i$dunnsidak_orthopvalues <- tmp_emp



j <- read.table("weigel_athaliana_IBE_final_orthogroup.txt", h = T)
species <- rep("athaliana", length(j$orthogroup))
j <- cbind(j,species)
j <- tested %>%
  left_join(j, by = c("orthogroup"))
j<-j[complete.cases(j), ]
pvals <- j$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
j$dunnsidak_orthopvalues <- tmp_emp


k <- read.table("weigel_capsella_final_orthogroup.txt", h = T)
species <- rep("capsella", length(k$orthogroup))
k <- cbind(k,species)
k <- tested %>%
  left_join(k, by = c("orthogroup"))
k<-k[complete.cases(k), ]
pvals <- k$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
k$dunnsidak_orthopvalues <- tmp_emp



l <- read.table("wright_atuber_final_orthogroup.txt", h = T)
species <- rep("atuber", length(l$orthogroup))
l <- cbind(l,species)
l <- tested %>%
  left_join(l, by = c("orthogroup"))
l<-l[complete.cases(l), ]
pvals <- l$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
l$dunnsidak_orthopvalues <- tmp_emp



m <- read.table("todesco_hann_final_orthogroup.txt", h = T)
species <- rep("hann", length(m$orthogroup))
m <- cbind(m,species)
m <- tested %>%
  left_join(m, by = c("orthogroup"))
m<-m[complete.cases(m), ]
pvals <- m$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
m$dunnsidak_orthopvalues <- tmp_emp



n <- read.table("todesco_hargo_final_orthogroup.txt", h = T)
species <- rep("hargo", length(n$orthogroup))
n <- cbind(n,species)
n <- tested %>%
  left_join(n, by = c("orthogroup"))
n<-n[complete.cases(n), ]
pvals <- n$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
n$dunnsidak_orthopvalues <- tmp_emp



o <- read.table("evans_ptricho_final_orthogroup.txt", h = T)
species <- rep("ptricho", length(o$orthogroup))
o <- cbind(o,species)
o <- tested %>%
  left_join(o, by = c("orthogroup"))
o<-o[complete.cases(o), ]
pvals <- o$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
o$dunnsidak_orthopvalues <- tmp_emp



p <- read.table("nocchi_bplaty_final_orthogroup.txt", h = T)
species <- rep("bplaty", length(p$orthogroup))
p <- cbind(p,species)
p <- tested %>%
  left_join(p, by = c("orthogroup"))
p<-p[complete.cases(p), ]
pvals <- p$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
p$dunnsidak_orthopvalues <- tmp_emp



q <- read.table("salojarvi_bpendula_final_orthogroup.txt", h = T)
species <- rep("bpendula", length(q$orthogroup))
q <- cbind(q,species)
q <- tested %>%
  left_join(q, by = c("orthogroup"))
q<-q[complete.cases(q), ]
pvals <- q$dunnsidak_orthopvalues
tmp_emp <- assign.pvalues(pvals)
q$dunnsidak_orthopvalues <- tmp_emp



results <- rbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q)
climate <- rep("mean_temp", length(results$species))

results <- cbind(results,climate)
colnames(results) <- c("Orthogroup","gene","min_CLR","scan_n","mean_emp_p","ortho_size", "ortho_DS", "species", "climate")
saveRDS(results, "orthogroup_results.rds")
