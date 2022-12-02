a <- read.table("kubota_ahalleri_final_orthogroup.txt", h = T)
species <- rep("ahalleri", length(a$orthogroup))
a <- cbind(a,species)

b <- read.table("ingvarsson_ptremula_final_orthogroup.txt", h = T)
species <- rep("ptremula", length(b$orthogroup))
b <- cbind(b,species)

c <- read.table("lowry_phalli_final_orthogroup.txt", h = T)
species <- rep("phalli", length(c$orthogroup))
c <- cbind(c,species)


d <- read.table("mitchell_bstricta_final_orthogroup.txt", h = T)
species <- rep("bstricta", length(d$orthogroup))
d <- cbind(d,species)

e <- read.table("murray_ealb_final_orthogroup.txt", h = T)
species <- rep("ealb", length(e$orthogroup))
e <- cbind(e,species)

f <- read.table("murray_esid_final_orthogroup.txt", h = T)
species <- rep("esid", length(f$orthogroup))
f <- cbind(f,species)

g <- read.table("murray_emag_final_orthogroup.txt", h = T)
species <- rep("emag", length(g$orthogroup))
g <- cbind(g,species)

h <- read.table("savolainen_alyrata_final_orthogroup.txt", h = T)
species <- rep("alyrata", length(h$orthogroup))
h <- cbind(h,species)


i <- read.table("tiffin_mtruncatula_final_orthogroup.txt", h = T)
species <- rep("mtruncatula", length(i$orthogroup))
i <- cbind(i,species)

j <- read.table("weigel_athaliana_IBE_final_orthogroup.txt", h = T)
species <- rep("athaliana", length(j$orthogroup))
j <- cbind(j,species)

k <- read.table("weiger_capsella_final_orthogroup.txt", h = T)
species <- rep("capsella", length(k$orthogroup))
k <- cbind(k,species)

l <- read.table("wright_atuber_final_orthogroup.txt", h = T)
species <- rep("atuber", length(l$orthogroup))
l <- cbind(l,species)


results <- rbind(a,b,c,d,e,f,g,h,i,j,k,l)
climate <- rep("mean_temp", length(results$species))

results <- cbind(results,climate)
colnames(results) <- c("gene","min_CLR","scan_n","mean_emp_p","Orthogroup","ortho_size", "ortho_DS", "species", "climate")
saveRDS(results, "orthogroup_results.rds")
