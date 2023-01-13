
setwd("/Users/gnocc/Desktop/PicMin/")
lib <- c("cowplot","poolr","mvmeta","qvalue","tidyr","ape","VGAM","ggExtra","pbmcapply","parallel","wCorr","data.table","ggplot2","viridis","ggridges","dplyr","readr")
sapply(lib,library,character.only=T)

# Function Library --------------------------------------------------------
source("picmin.R")


############################################################################
# Assemble the Orthofinder Outputs to Analyse
n_cores = 1 # Number of cpu
orthogroup_cutoff <- 7 # Minimum number of species in an OG
# max_paralog_per_OG <- 10

# Fetch the processed Orthogroup Pvals
OG_pvals <- readRDS("orthogroup_results.rds")
OG_pvals <- as.data.table(OG_pvals)

##################################################

# Define climate variaibles here
climate_vars <- unique(OG_pvals$climate)
focal_climate <- "mean_temp"
OG_maxP <- OG_pvals[climate == focal_climate,]
focal_datasets_climate <- unique(OG_pvals$species)



# Remove orthogroups with not enough datasets...
remove_small_OG = names(table(OG_maxP$Orthogroup)[table(OG_maxP$Orthogroup) < orthogroup_cutoff])

# Remove orthogroups with too many paralogs in any species...
# remove_paralogs = unique(OG_maxP[Ngenes_per_species > max_paralog_per_OG,Orthogroup])

# Filter
OG_maxP_test <- OG_maxP[!OG_maxP$Orthogroup %in% unique(remove_small_OG),]

# And get minP
minPval_OG <- OG_maxP_test[,.(minP=min(ortho_DS)),by=Orthogroup]

picmin_cor_mats <- pbmclapply(3:length(focal_datasets_climate),function(n){

  ee <- array (NA,c(10000,n))
  for (i in 1:n){
    temp <- runif(10000)
    ee[,i] <- EmpiricalPs(temp)
  }
  emp_p_null_dat30_unif <- ee

  for (kk in 1:nrow(emp_p_null_dat30_unif)){
    ind1 <- order(emp_p_null_dat30_unif[kk,])
    emp_p_null_dat30_unif [kk,] <- emp_p_null_dat30_unif[kk,ind1]
  }

  cor_ord_unif <- cor (emp_p_null_dat30_unif[,2:n])
  cor_ord_unif
},mc.cores=n_cores)


names(picmin_cor_mats) <- as.character(3:length(focal_datasets_climate))


picmin_res <- rbindlist(pbmclapply(unique(OG_maxP_test$Orthogroup),function(OG){
  set.seed(1000)
  tmp <- OG_maxP_test[OG_maxP_test$Orthogroup == OG,]
  # Run first with dunn-sidak correction
  # first_run = PicMin_bugfix(tmp$max_sdP_DS, picmin_cor_mats[[as.character(nrow(tmp))]], numReps = 100,correction = "dunn-sidak")
  first_run = PicMin_bugfix(tmp$ortho_DS, picmin_cor_mats[[as.character(nrow(tmp))]], numReps = 1000,correction = "tippett")

  if(first_run$p > 0.1){
    return(first_run)
  } else {
    # # How many iterations to run Tippett
    second_run_iter = ifelse(10/first_run$p > 1e8, 1e8, round(10/first_run$p)*10)
    # # Run
    second_run = PicMin_bugfix(tmp$ortho_DS, picmin_cor_mats[[as.character(nrow(tmp))]], numReps = second_run_iter,correction = "tippett")
    return(second_run)
  }
},mc.cores=n_cores))


picmin_res$picmin_fdr <- p.adjust(picmin_res$p,method = "fdr")
picmin_res$Orthogroup <- unique(OG_maxP_test$Orthogroup)
picmin_res$climate_var <- focal_climate



  picmin_res[picmin_res$p == min(picmin_res$p),]
  picmin_res[order(p),]

  # Merge with the minP as well
  picmin_res = merge(picmin_res,minPval_OG,by = "Orthogroup")

  # Save these
  picmin_res_output <- list(OG_maxP=OG_maxP_test,
                            picmin_res=picmin_res)
  saveRDS(picmin_res_output, "gab_picmin_results.rds")
