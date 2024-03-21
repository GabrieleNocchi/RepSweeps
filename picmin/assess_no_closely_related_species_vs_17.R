a <- readRDS("gab_picmin_results_1.rds")
a <- a$picmin_res
a <- a[a$picmin_fdr < 0.5,]
a <- a$Orthogroup



b <- readRDS("gab_picmin_results_2.rds")
b <- b$picmin_res
b <- b[b$picmin_fdr < 0.5,]
b <- b$Orthogroup


c <- readRDS("gab_picmin_results_3.rds")
c <- c$picmin_res
c <- c[c$picmin_fdr < 0.5,]
c <- c$Orthogroup


d <- readRDS("gab_picmin_results_4.rds")
d <- d$picmin_res
d <- d[d$picmin_fdr < 0.5,]
d <- d$Orthogroup



e <- readRDS("gab_picmin_results_5.rds")
e <- e$picmin_res
e <- e[e$picmin_fdr < 0.5,]
e <- e$Orthogroup


f <- readRDS("gab_picmin_results_6.rds")
f <- f$picmin_res
f <- f[f$picmin_fdr < 0.5,]
f <- f$Orthogroup


g <- readRDS("gab_picmin_results_7.rds")
g <- g$picmin_res
g <- g[g$picmin_fdr < 0.5,]
g <- g$Orthogroup


h <- readRDS("gab_picmin_results_8.rds")
h <- h$picmin_res
h <- h[h$picmin_fdr < 0.5,]
h <- h$Orthogroup


i <- readRDS("gab_picmin_results_9.rds")
i <- i$picmin_res
i <- i[i$picmin_fdr < 0.5,]
i <- i$Orthogroup


separated <- c(a,b,c,d,e,f,g,h,i)
separated <- unique(sort(separated))
separated <- as.data.frame(separated)
colnames(separated) <- "Orthogroup"

tot <- readRDS("gab_picmin_results.rds")
tot <- tot$picmin_res
tot <- tot[tot$picmin_fdr < 0.5,]



library(dplyr)
new <- separated %>%
  left_join(tot, by = c("Orthogroup"))
  new <-new[complete.cases(new), ]

saveRDS(new, "gab_picmin_results_overlap_17_0.5_13_0.5.rds")
