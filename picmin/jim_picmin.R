# Orthogroup-level convergence using PicMin
# This is a variant where alpha_a is a flexible parameter and calculated according to the top 1000 orthogroups...
lib <- c("cowplot","poolr","mvmeta","qvalue","tidyr","ape","VGAM","ggExtra","pbmcapply","parallel","wCorr","data.table","ggplot2","viridis","ggridges","dplyr","readr")
sapply(lib,library,character.only=T)

# Function Library --------------------------------------------------------
source("R/PicMin.R")

############################################################################
# What GEA results are we looking at?
run_name = "220927"
# picmin_reps = 1000000
output_name = "COMBINED25species_NOMIRROR_221026"
pvals_file = paste0("outputs/GEA_res/run",run_name,"_",output_name,"_WZA_OG_pvals.rds")

############################################################################
# Assemble the Orthofinder Outputs to Analyse
n_cores = 5 # Number of cpu
orthogroup_cutoff <- 20 # Minimum number of species in an OG
# max_paralog_per_OG <- 10

# Fetch the processed Orthogroup Pvals
OG_pvals <- readRDS(pvals_file)

############################################################################################################################################################
# Prepare GEA data -------------------------------------------------------
focal_datasets <- list.files("outputs/GEA_res",pattern = run_name)
focal_datasets <- grep("processed",focal_datasets,invert = T,value = T)
focal_datasets <- grep("CoAdapTree",focal_datasets,invert = T,value = T)

# Define climate variaibles here
climate_vars <- unique(OG_pvals$climate)

# We loop over our focal climates to calculate orthogroup-level results for all climates
for(focal_climate in climate_vars){

  print(paste0(">>> STARTING PICMIN FOR:", focal_climate))

  # Subset the OG for climate res of focus
  OG_maxP <- OG_pvals[climate == focal_climate,]
  focal_datasets_climate <- unique(OG_pvals$species)

  # Remove orthogroups with not enough datasets...
  remove_small_OG = names(table(OG_maxP$Orthogroup)[table(OG_maxP$Orthogroup) < orthogroup_cutoff])

  # Remove orthogroups with too many paralogs in any species...
  # remove_paralogs = unique(OG_maxP[Ngenes_per_species > max_paralog_per_OG,Orthogroup])

  # Filter
  OG_maxP_test <- OG_maxP[!OG_maxP$Orthogroup %in% unique(remove_small_OG),]

  # And get minP
  minPval_OG <- OG_maxP_test[,.(minP=min(min_sdP_DS)),by=Orthogroup]

  # Make correlation matrices for Tippett's of size 3-N
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


  # Run PicMin - Time testing -----------------------------------------------
  message(paste0(">>> Starting PicMin/FixMin"))

  picmin_res <- rbindlist(pbmclapply(unique(OG_maxP_test$Orthogroup),function(OG){
    set.seed(1000)
    tmp <- OG_maxP_test[OG_maxP_test$Orthogroup == OG,]
    # Run first with dunn-sidak correction
    # first_run = PicMin_bugfix(tmp$max_sdP_DS, picmin_cor_mats[[as.character(nrow(tmp))]], numReps = 100,correction = "dunn-sidak")
    first_run = PicMin_bugfix(tmp$min_sdP_DS, picmin_cor_mats[[as.character(nrow(tmp))]], numReps = 1000,correction = "tippett")

    if(first_run$p > 0.1){
      return(first_run)
    } else {
      # # How many iterations to run Tippett
      second_run_iter = ifelse(10/first_run$p > 1e8, 1e8, round(10/first_run$p)*10)
      # # Run
      second_run = PicMin_bugfix(tmp$min_sdP_DS, picmin_cor_mats[[as.character(nrow(tmp))]], numReps = second_run_iter,correction = "tippett")
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
  saveRDS(picmin_res_output,
          paste0("outputs/picmin_res_allgenes_minSpecies",orthogroup_cutoff,"_run",run_name,"_",focal_climate,"_sdpvals_",output_name,".rds"))

}

# Quick check outputs
for(focal_climate in climate_vars){
  tmp = readRDS(paste0("outputs/picmin_res_allgenes_minSpecies",orthogroup_cutoff,"_run",run_name,"_",focal_climate,"_sdpvals_",output_name,".rds"))
  N_outliers = nrow(tmp$picmin_res[picmin_fdr < 0.5,])
  print(paste0(focal_climate,": ",N_outliers))
}
#
# earlier_run = readRDS("outputs/test_result_sets_for_different_picmin_runs_annual_precip.rds")
#
# merge_res = merge(test,earlier_run$fixmin[,.(Orthogroup,p)],by="Orthogroup")
#
# ggplot(merge_res,aes(-log10(p.adjust(p.x,"fdr")),-log10(p.adjust(p.y,"fdr")))) + geom_point() +
#   geom_abline(colour = "red2")+
#   stat_density2d_filled(alpha = 0.5)
