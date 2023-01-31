library(ComplexHeatmap)

assign.rank <- function(array){
  pvalues <- array(0, length(array))

  ordered.indexes <- order(array, decreasing = TRUE)

  j <- length(array)
  for( i in ordered.indexes ){
    pvalues[i] <- j/length(array)
    j <- j-1
  }

  return(pvalues)
}
##### change threshold accordingly
thr <- 0.01

##### Preparing Gabriele top thr % orthogroup p values
##### wright_atuber
gab <- read.table("wright_atuber_final_orthogroup.txt", header=TRUE)
gab_ortho_p <- gab[,7]
gab_ranks <- assign.rank(gab_ortho_p)
gab <- cbind(gab,gab_ranks)
gab_outliers <- gab[gab$gab_ranks < thr,]

##### Preparing Jim top thr % ortho p values
jim <- readRDS("run220927_COMBINED25species_updatedOF_221213_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Amaranthus tuberculatus",]
jim <-jim[order(jim$Orthogroup),]
jim_ranks <- sapply(split(jim$min_sdP_DS, jim$climate), assign.rank)
jim_ranks <- data.frame(jim_ranks)

 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 annual_precip <- annual_precip[annual_precip$X2 < thr,]

 isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 isothermality <- isothermality[isothermality$X2 < thr,]

 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X2 < thr,]

 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_diurnal <- mean_diurnal[mean_diurnal$X2 < thr,]

 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp <- mean_temp[mean_temp$X2 < thr,]

 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X2 < thr,]

 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X2 < thr,]

 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X2 < thr,]

 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X2 < thr,]

 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X2 < thr,]

 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 prec_clim_change <- prec_clim_change[prec_clim_change$X2 < thr,]

 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X2 < thr,]

 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_month <- precip_dry_month[precip_dry_month$X2 < thr,]

 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X2 < thr,]

 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_seasonality <- precip_seasonality[precip_seasonality$X2 < thr,]

 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X2 < thr,]

 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_month <- precip_wet_month[precip_wet_month$X2 < thr,]

 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X2 < thr,]

 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_range <- temp_range[temp_range$X2 < thr,]

 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 temp_seasonality <- temp_seasonality[temp_seasonality$X2 < thr,]

 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))
 tmax_clim_change <- tmax_clim_change[tmax_clim_change$X2 < thr,]

 jim_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)


  res <- lapply(jim_list, function(x){
  test <- length(intersect(gab_outliers$orthogroup, x$X1))
  return(test)
 })

 # Hypergeometric expectancy
 #res_percent <- lapply(res,"/",length(annual_precip$X1))
 total <- length(unique(sort(jim$Orthogroup)))
 percent <- 100/(thr * 100)
 expect_outliers <- total/percent
 hyper_exp <- (expect_outliers * expect_outliers)/total
 res_percent <- lapply(res,"/",hyper_exp)

df <- data.frame(res_percent)
colnames(df) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df) <- c("wright_atuber")






##### Preparing Gabriele top thr % orthogroup p values
###### lowry_phalli
gab <- read.table("lowry_phallii_final_orthogroup.txt", header=TRUE)
gab_ortho_p <- gab[,7]
gab_ranks <- assign.rank(gab_ortho_p)
gab <- cbind(gab,gab_ranks)
gab_outliers <- gab[gab$gab_ranks < thr,]

###### Preparing Jim top thr % ortho p values
jim <- readRDS("run220927_COMBINED25species_updatedOF_221213_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Panicum hallii",]
jim <-jim[order(jim$Orthogroup),]
jim_ranks <- sapply(split(jim$min_sdP_DS, jim$climate), assign.rank)
jim_ranks <- data.frame(jim_ranks)

 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 annual_precip <- annual_precip[annual_precip$X2 < thr,]

 isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 isothermality <- isothermality[isothermality$X2 < thr,]

 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X2 < thr,]

 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_diurnal <- mean_diurnal[mean_diurnal$X2 < thr,]

 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp <- mean_temp[mean_temp$X2 < thr,]

 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X2 < thr,]

 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X2 < thr,]

 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X2 < thr,]

 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X2 < thr,]

 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X2 < thr,]

 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 prec_clim_change <- prec_clim_change[prec_clim_change$X2 < thr,]

 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X2 < thr,]

 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_month <- precip_dry_month[precip_dry_month$X2 < thr,]

 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X2 < thr,]

 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_seasonality <- precip_seasonality[precip_seasonality$X2 < thr,]

 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X2 < thr,]

 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_month <- precip_wet_month[precip_wet_month$X2 < thr,]

 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X2 < thr,]

 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_range <- temp_range[temp_range$X2 < thr,]

 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 temp_seasonality <- temp_seasonality[temp_seasonality$X2 < thr,]

 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))
 tmax_clim_change <- tmax_clim_change[tmax_clim_change$X2 < thr,]

 jim_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)


 res <- lapply(jim_list, function(x){
  test <- length(intersect(gab_outliers$orthogroup, x$X1))
  return(test)
 })

 # Hypergeometric expectancy
 #res_percent <- lapply(res,"/",length(annual_precip$X1))
total <- length(unique(sort(jim$Orthogroup)))
 percent <- 100/(thr * 100)
 expect_outliers <- total/percent
 hyper_exp <- (expect_outliers * expect_outliers)/total
 res_percent <- lapply(res,"/",hyper_exp)

df2 <- data.frame(res_percent)
colnames(df2) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df2) <- c("lowry_phallii")






 ##### Preparing Gabriele top thr % orthogroup p values
 ##### weigel_athaliana
 gab <- read.table("weigel_athaliana_IBE_final_orthogroup.txt", header=TRUE)
 gab_ortho_p <- gab[,7]
 gab_ranks <- assign.rank(gab_ortho_p)
 gab <- cbind(gab,gab_ranks)
 gab_outliers <- gab[gab$gab_ranks < thr,]

##### Preparing Jim top thr % ortho p values
 jim <- readRDS("run220927_COMBINED25species_updatedOF_221213_WZA_OG_pvals.rds")
 jim <- jim[jim$species == "Arabidopsis thaliana",]
 jim <-jim[order(jim$Orthogroup),]
 jim_ranks <- sapply(split(jim$min_sdP_DS, jim$climate), assign.rank)
 jim_ranks <- data.frame(jim_ranks)

  annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
  annual_precip <- annual_precip[annual_precip$X2 < thr,]

  isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
  isothermality <- isothermality[isothermality$X2 < thr,]

  max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
  max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X2 < thr,]

  mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
  mean_diurnal <- mean_diurnal[mean_diurnal$X2 < thr,]

  mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
  mean_temp <- mean_temp[mean_temp$X2 < thr,]

  mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
  mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X2 < thr,]

  mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
  mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X2 < thr,]

  mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
  mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X2 < thr,]

  mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
  mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X2 < thr,]

  min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
  min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X2 < thr,]

  prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
  prec_clim_change <- prec_clim_change[prec_clim_change$X2 < thr,]

  precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
  precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X2 < thr,]

  precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
  precip_dry_month <- precip_dry_month[precip_dry_month$X2 < thr,]

  precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
  precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X2 < thr,]

  precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
  precip_seasonality <- precip_seasonality[precip_seasonality$X2 < thr,]

  precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
  precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X2 < thr,]

  precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
  precip_wet_month <- precip_wet_month[precip_wet_month$X2 < thr,]

  precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
  precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X2 < thr,]

  temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
  temp_range <- temp_range[temp_range$X2 < thr,]

  temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
  temp_seasonality <- temp_seasonality[temp_seasonality$X2 < thr,]

  tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))
  tmax_clim_change <- tmax_clim_change[tmax_clim_change$X2 < thr,]

  jim_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)


  res <- lapply(jim_list, function(x){
   test <- length(intersect(gab_outliers$orthogroup, x$X1))
   return(test)
  })

  # Hypergeometric expectancy
  #res_percent <- lapply(res,"/",length(annual_precip$X1))
  total <- length(unique(sort(jim$Orthogroup)))
  percent <- 100/(thr * 100)
  expect_outliers <- total/percent
  hyper_exp <- (expect_outliers * expect_outliers)/total
  res_percent <- lapply(res,"/",hyper_exp)

 df3 <- data.frame(res_percent)
 colnames(df3) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
 rownames(df3) <- c("weigel_athaliana")






##### Preparing Gabriele top thr % orthogroup p values
##### weiger_capsella
gab <- read.table("weigel_capsella_final_orthogroup.txt", header=TRUE)
gab_ortho_p <- gab[,7]
gab_ranks <- assign.rank(gab_ortho_p)
gab <- cbind(gab,gab_ranks)
gab_outliers <- gab[gab$gab_ranks < thr,]

##### Preparing Jim top thr % ortho p values
jim <- readRDS("run220927_COMBINED25species_updatedOF_221213_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Capsella rubella",]
jim <-jim[order(jim$Orthogroup),]
jim_ranks <- sapply(split(jim$min_sdP_DS, jim$climate), assign.rank)
jim_ranks <- data.frame(jim_ranks)

 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 annual_precip <- annual_precip[annual_precip$X2 < thr,]

 isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 isothermality <- isothermality[isothermality$X2 < thr,]

 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X2 < thr,]

 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_diurnal <- mean_diurnal[mean_diurnal$X2 < thr,]

 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp <- mean_temp[mean_temp$X2 < thr,]

 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X2 < thr,]

 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X2 < thr,]

 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X2 < thr,]

 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X2 < thr,]

 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X2 < thr,]

 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 prec_clim_change <- prec_clim_change[prec_clim_change$X2 < thr,]

 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X2 < thr,]

 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_month <- precip_dry_month[precip_dry_month$X2 < thr,]

 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X2 < thr,]

 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_seasonality <- precip_seasonality[precip_seasonality$X2 < thr,]

 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X2 < thr,]

 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_month <- precip_wet_month[precip_wet_month$X2 < thr,]

 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X2 < thr,]

 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_range <- temp_range[temp_range$X2 < thr,]

 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 temp_seasonality <- temp_seasonality[temp_seasonality$X2 < thr,]

 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))
 tmax_clim_change <- tmax_clim_change[tmax_clim_change$X2 < thr,]

 jim_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)


 res <- lapply(jim_list, function(x){
  test <- length(intersect(gab_outliers$orthogroup, x$X1))
  return(test)
 })

 # Hypergeometric expectancy
 #res_percent <- lapply(res,"/",length(annual_precip$X1))
total <- length(unique(sort(jim$Orthogroup)))
 percent <- 100/(thr * 100)
 expect_outliers <- total/percent
 hyper_exp <- (expect_outliers * expect_outliers)/total
 res_percent <- lapply(res,"/",hyper_exp)

df4 <- data.frame(res_percent)
colnames(df4) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df4) <- c("weigel_capsella")






##### Preparing Gabriele top thr % orthogroup p values
###### murray_emag
gab <- read.table("murray_emag_final_orthogroup.txt", header=TRUE)
gab_ortho_p <- gab[,7]
gab_ranks <- assign.rank(gab_ortho_p)
gab <- cbind(gab,gab_ranks)
gab_outliers <- gab[gab$gab_ranks < thr,]

##### Preparing Jim top thr % ortho p values
jim <- readRDS("run220927_COMBINED25species_updatedOF_221213_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Eucalyptus magnificata",]
jim <-jim[order(jim$Orthogroup),]
jim_ranks <- sapply(split(jim$min_sdP_DS, jim$climate), assign.rank)
jim_ranks <- data.frame(jim_ranks)

 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 annual_precip <- annual_precip[annual_precip$X2 < thr,]

 #isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 #isothermality <- isothermality[isothermality$X2 < thr,]
 isothermality <- data.frame(X1=c(rep("gabriele",1255)), X2 = c(rep(1,1255)))

 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X2 < thr,]

 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_diurnal <- mean_diurnal[mean_diurnal$X2 < thr,]

 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp <- mean_temp[mean_temp$X2 < thr,]

 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X2 < thr,]

 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X2 < thr,]

 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X2 < thr,]

 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X2 < thr,]

 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X2 < thr,]

 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 prec_clim_change <- prec_clim_change[prec_clim_change$X2 < thr,]

 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X2 < thr,]

 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_month <- precip_dry_month[precip_dry_month$X2 < thr,]

 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X2 < thr,]

 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_seasonality <- precip_seasonality[precip_seasonality$X2 < thr,]

 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X2 < thr,]

 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_month <- precip_wet_month[precip_wet_month$X2 < thr,]

 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X2 < thr,]

 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_range <- temp_range[temp_range$X2 < thr,]

 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 temp_seasonality <- temp_seasonality[temp_seasonality$X2 < thr,]

 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))
 tmax_clim_change <- tmax_clim_change[tmax_clim_change$X2 < thr,]

 jim_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)


 res <- lapply(jim_list, function(x){
  test <- length(intersect(gab_outliers$orthogroup, x$X1))
  return(test)
 })

 # Hypergeometric expectancy
 #res_percent <- lapply(res,"/",length(annual_precip$X1))
 total <- length(unique(sort(jim$Orthogroup)))
 percent <- 100/(thr * 100)
 expect_outliers <- total/percent
 hyper_exp <- (expect_outliers * expect_outliers)/total
 res_percent <- lapply(res,"/",hyper_exp)

df5 <- data.frame(res_percent)
colnames(df5) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df5) <- c("murray_emag")






##### Preparing Gabriele top thr % orthogroup p values
###### kubota_ahalleri
gab <- read.table("kubota_ahalleri_final_orthogroup.txt", header=TRUE)
gab_ortho_p <- gab[,7]
gab_ranks <- assign.rank(gab_ortho_p)
gab <- cbind(gab,gab_ranks)
gab_outliers <- gab[gab$gab_ranks < thr,]

##### Preparing Jim top thr % ortho p values
jim <- readRDS("run220927_COMBINED25species_updatedOF_221213_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Arabidopsis halleri",]
jim <-jim[order(jim$Orthogroup),]
jim_ranks <- sapply(split(jim$min_sdP_DS, jim$climate), assign.rank)
jim_ranks <- data.frame(jim_ranks)

 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 annual_precip <- annual_precip[annual_precip$X2 < thr,]

 isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 isothermality <- isothermality[isothermality$X2 < thr,]

 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X2 < thr,]

 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_diurnal <- mean_diurnal[mean_diurnal$X2 < thr,]

 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp <- mean_temp[mean_temp$X2 < thr,]

 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X2 < thr,]

 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X2 < thr,]

 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X2 < thr,]

 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X2 < thr,]

 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X2 < thr,]

 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 prec_clim_change <- prec_clim_change[prec_clim_change$X2 < thr,]

 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X2 < thr,]

 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_month <- precip_dry_month[precip_dry_month$X2 < thr,]

 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X2 < thr,]

 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_seasonality <- precip_seasonality[precip_seasonality$X2 < thr,]

 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X2 < thr,]

 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_month <- precip_wet_month[precip_wet_month$X2 < thr,]

 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X2 < thr,]

 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_range <- temp_range[temp_range$X2 < thr,]

 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 temp_seasonality <- temp_seasonality[temp_seasonality$X2 < thr,]

 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))
 tmax_clim_change <- tmax_clim_change[tmax_clim_change$X2 < thr,]

 jim_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)


 res <- lapply(jim_list, function(x){
  test <- length(intersect(gab_outliers$orthogroup, x$X1))
  return(test)
 })

 # Hypergeometric expectancy
 #res_percent <- lapply(res,"/",length(annual_precip$X1))
 total <- length(unique(sort(jim$Orthogroup)))
 percent <- 100/(thr * 100)
 expect_outliers <- total/percent
 hyper_exp <- (expect_outliers * expect_outliers)/total
 res_percent <- lapply(res,"/",hyper_exp)

df6 <- data.frame(res_percent)
colnames(df6) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df6) <- c("kubota_ahalleri")






##### Preparing Gabriele top thr % orthogroup p values
###### murray_esid
gab <- read.table("murray_esid_final_orthogroup.txt", header=TRUE)
gab_ortho_p <- gab[,7]
gab_ranks <- assign.rank(gab_ortho_p)
gab <- cbind(gab,gab_ranks)
gab_outliers <- gab[gab$gab_ranks < thr,]

##### Preparing Jim top thr % ortho p values
jim <- readRDS("run220927_COMBINED25species_updatedOF_221213_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Eucalyptus sideroxylon",]
jim <-jim[order(jim$Orthogroup),]
jim_ranks <- sapply(split(jim$min_sdP_DS, jim$climate), assign.rank)
jim_ranks <- data.frame(jim_ranks)

 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 annual_precip <- annual_precip[annual_precip$X2 < thr,]

 isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 isothermality <- isothermality[isothermality$X2 < thr,]

 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X2 < thr,]

 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_diurnal <- mean_diurnal[mean_diurnal$X2 < thr,]

 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp <- mean_temp[mean_temp$X2 < thr,]

 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X2 < thr,]

 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X2 < thr,]

 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X2 < thr,]

 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X2 < thr,]

 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X2 < thr,]

 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 prec_clim_change <- prec_clim_change[prec_clim_change$X2 < thr,]

 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X2 < thr,]

 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_month <- precip_dry_month[precip_dry_month$X2 < thr,]

 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X2 < thr,]

 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_seasonality <- precip_seasonality[precip_seasonality$X2 < thr,]

 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X2 < thr,]

 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_month <- precip_wet_month[precip_wet_month$X2 < thr,]

 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X2 < thr,]

 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_range <- temp_range[temp_range$X2 < thr,]

 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 temp_seasonality <- temp_seasonality[temp_seasonality$X2 < thr,]

 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))
 tmax_clim_change <- tmax_clim_change[tmax_clim_change$X2 < thr,]

 jim_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)


 res <- lapply(jim_list, function(x){
  test <- length(intersect(gab_outliers$orthogroup, x$X1))
  return(test)
 })

 # Hypergeometric expectancy
 #res_percent <- lapply(res,"/",length(annual_precip$X1))
 total <- length(unique(sort(jim$Orthogroup)))
 percent <- 100/(thr * 100)
 expect_outliers <- total/percent
 hyper_exp <- (expect_outliers * expect_outliers)/total
 res_percent <- lapply(res,"/",hyper_exp)

df7 <- data.frame(res_percent)
colnames(df7) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df7) <- c("murray_esid")






# ##### Preparing Gabriele top thr % orthogroup p values
# ###### savolainen_alyrata
# gab <- read.table("savolainen_alyrata_final_orthogroup.txt", header=TRUE)
# gab_ortho_p <- gab[,7]
# gab_ranks <- assign.rank(gab_ortho_p)
# gab <- cbind(gab,gab_ranks)
# gab_outliers <- gab[gab$gab_ranks < thr,]
#
# ##### Preparing Jim top thr % ortho p values
# jim <- readRDS("run220927_COMBINED25species_updatedOF_221213_WZA_OG_pvals.rds")
# jim <- jim[jim$species == "Arabidopsis lyrata",]
# jim <-jim[order(jim$Orthogroup),]
# jim_ranks <- sapply(split(jim$min_sdP_DS, jim$climate), assign.rank)
# jim_ranks <- data.frame(jim_ranks)
#
#  annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
#  annual_precip <- annual_precip[annual_precip$X2 < thr,]
#
#  isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
#  isothermality <- isothermality[isothermality$X2 < thr,]
#
#  max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
#  max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X2 < thr,]
#
#  mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
#  mean_diurnal <- mean_diurnal[mean_diurnal$X2 < thr,]
#
#  mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
#  mean_temp <- mean_temp[mean_temp$X2 < thr,]
#
#  mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
#  mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X2 < thr,]
#
#  mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
#  mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X2 < thr,]
#
#  mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
#  mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X2 < thr,]
#
#  mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
#  mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X2 < thr,]
#
#  min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
#  min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X2 < thr,]
#
#  prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
#  prec_clim_change <- prec_clim_change[prec_clim_change$X2 < thr,]
#
#  precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
#  precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X2 < thr,]
#
#  precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
#  precip_dry_month <- precip_dry_month[precip_dry_month$X2 < thr,]
#
#  precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
#  precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X2 < thr,]
#
#  precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
#  precip_seasonality <- precip_seasonality[precip_seasonality$X2 < thr,]
#
#  precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
#  precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X2 < thr,]
#
#  precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
#  precip_wet_month <- precip_wet_month[precip_wet_month$X2 < thr,]
#
#  precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
#  precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X2 < thr,]
#
#  temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
#  temp_range <- temp_range[temp_range$X2 < thr,]
#
#  temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
#  temp_seasonality <- temp_seasonality[temp_seasonality$X2 < thr,]
#
#  tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))
#  tmax_clim_change <- tmax_clim_change[tmax_clim_change$X2 < thr,]
#
#  jim_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)
#
#
#  res <- lapply(jim_list, function(x){
#   test <- length(intersect(gab_outliers$orthogroup, x$X1))
#   return(test)
#  })
#
#  # Hypergeometric expectancy
#  #res_percent <- lapply(res,"/",length(annual_precip$X1))
#  total <- length(unique(sort(jim$Orthogroup)))
#  percent <- 100/(thr * 100)
#  expect_outliers <- total/percent
#  hyper_exp <- (expect_outliers * expect_outliers)/total
#  res_percent <- lapply(res,"/",hyper_exp)
#
# df8 <- data.frame(res_percent)
# colnames(df8) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
# rownames(df8) <- c("savolainen_alyrata")






##### Preparing Gabriele top thr % orthogroup p values
##### mitchell_bstricta
gab <- read.table("mitchell_bstricta_final_orthogroup.txt", header=TRUE)
gab_ortho_p <- gab[,7]
gab_ranks <- assign.rank(gab_ortho_p)
gab <- cbind(gab,gab_ranks)
gab_outliers <- gab[gab$gab_ranks < thr,]

##### Preparing Jim top thr % ortho p values
jim <- readRDS("run220927_COMBINED25species_updatedOF_221213_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Boechera stricta",]
jim <-jim[order(jim$Orthogroup),]
jim_ranks <- sapply(split(jim$min_sdP_DS, jim$climate), assign.rank)
jim_ranks <- data.frame(jim_ranks)

 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 annual_precip <- annual_precip[annual_precip$X2 < thr,]

 isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 isothermality <- isothermality[isothermality$X2 < thr,]

 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X2 < thr,]

 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_diurnal <- mean_diurnal[mean_diurnal$X2 < thr,]

 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp <- mean_temp[mean_temp$X2 < thr,]

 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X2 < thr,]

 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X2 < thr,]

 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X2 < thr,]

 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X2 < thr,]

 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X2 < thr,]

 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 prec_clim_change <- prec_clim_change[prec_clim_change$X2 < thr,]

 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X2 < thr,]

 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_month <- precip_dry_month[precip_dry_month$X2 < thr,]

 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X2 < thr,]

 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_seasonality <- precip_seasonality[precip_seasonality$X2 < thr,]

 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X2 < thr,]

 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_month <- precip_wet_month[precip_wet_month$X2 < thr,]

 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X2 < thr,]

 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_range <- temp_range[temp_range$X2 < thr,]

 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 temp_seasonality <- temp_seasonality[temp_seasonality$X2 < thr,]

 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))
 tmax_clim_change <- tmax_clim_change[tmax_clim_change$X2 < thr,]

 jim_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)


 res <- lapply(jim_list, function(x){
  test <- length(intersect(gab_outliers$orthogroup, x$X1))
  return(test)
 })

 # Hypergeometric expectancy
 #res_percent <- lapply(res,"/",length(annual_precip$X1))
 total <- length(unique(sort(jim$Orthogroup)))
 percent <- 100/(thr * 100)
 expect_outliers <- total/percent
 hyper_exp <- (expect_outliers * expect_outliers)/total
 res_percent <- lapply(res,"/",hyper_exp)

df9 <- data.frame(res_percent)
colnames(df9) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df9) <- c("mitchell_bstricta")






##### Preparing Gabriele top thr % orthogroup p values
##### murray_ealb
gab <- read.table("murray_ealb_final_orthogroup.txt", header=TRUE)
gab_ortho_p <- gab[,7]
gab_ranks <- assign.rank(gab_ortho_p)
gab <- cbind(gab,gab_ranks)
gab_outliers <- gab[gab$gab_ranks < thr,]

##### Preparing Jim top thr % ortho p values
jim <- readRDS("run220927_COMBINED25species_updatedOF_221213_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Eucalyptus albens",]
jim <-jim[order(jim$Orthogroup),]
jim_ranks <- sapply(split(jim$min_sdP_DS, jim$climate), assign.rank)
jim_ranks <- data.frame(jim_ranks)

 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 annual_precip <- annual_precip[annual_precip$X2 < thr,]

 isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 isothermality <- isothermality[isothermality$X2 < thr,]

 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X2 < thr,]

 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_diurnal <- mean_diurnal[mean_diurnal$X2 < thr,]

 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp <- mean_temp[mean_temp$X2 < thr,]

 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X2 < thr,]

 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X2 < thr,]

 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X2 < thr,]

 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X2 < thr,]

 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X2 < thr,]

 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 prec_clim_change <- prec_clim_change[prec_clim_change$X2 < thr,]

 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X2 < thr,]

 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_month <- precip_dry_month[precip_dry_month$X2 < thr,]

 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X2 < thr,]

 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_seasonality <- precip_seasonality[precip_seasonality$X2 < thr,]

 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X2 < thr,]

 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_month <- precip_wet_month[precip_wet_month$X2 < thr,]

 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X2 < thr,]

 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_range <- temp_range[temp_range$X2 < thr,]

 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 temp_seasonality <- temp_seasonality[temp_seasonality$X2 < thr,]

 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))
 tmax_clim_change <- tmax_clim_change[tmax_clim_change$X2 < thr,]

 jim_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)


 res <- lapply(jim_list, function(x){
  test <- length(intersect(gab_outliers$orthogroup, x$X1))
  return(test)
 })

 # Hypergeometric expectancy
 #res_percent <- lapply(res,"/",length(annual_precip$X1))
 total <- length(unique(sort(jim$Orthogroup)))
 percent <- 100/(thr * 100)
 expect_outliers <- total/percent
 hyper_exp <- (expect_outliers * expect_outliers)/total
 res_percent <- lapply(res,"/",hyper_exp)

df10 <- data.frame(res_percent)
colnames(df10) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df10) <- c("murray_ealb")






##### Preparing Gabriele top thr % orthogroup p values
###### tiffin_truncatula
gab <- read.table("tiffin_mtruncatula_final_orthogroup.txt", header=TRUE)
gab_ortho_p <- gab[,7]
gab_ranks <- assign.rank(gab_ortho_p)
gab <- cbind(gab,gab_ranks)
gab_outliers <- gab[gab$gab_ranks < thr,]

##### Preparing Jim top thr % ortho p values
jim <- readRDS("run220927_COMBINED25species_updatedOF_221213_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Medicago truncatula",]
jim <-jim[order(jim$Orthogroup),]
jim_ranks <- sapply(split(jim$min_sdP_DS, jim$climate), assign.rank)
jim_ranks <- data.frame(jim_ranks)

 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 annual_precip <- annual_precip[annual_precip$X2 < thr,]

 isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 isothermality <- isothermality[isothermality$X2 < thr,]

 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X2 < thr,]

 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_diurnal <- mean_diurnal[mean_diurnal$X2 < thr,]

 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp <- mean_temp[mean_temp$X2 < thr,]

 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X2 < thr,]

 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X2 < thr,]

 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X2 < thr,]

 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X2 < thr,]

 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X2 < thr,]

 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 prec_clim_change <- prec_clim_change[prec_clim_change$X2 < thr,]

 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X2 < thr,]

 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_month <- precip_dry_month[precip_dry_month$X2 < thr,]

 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X2 < thr,]

 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_seasonality <- precip_seasonality[precip_seasonality$X2 < thr,]

 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X2 < thr,]

 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_month <- precip_wet_month[precip_wet_month$X2 < thr,]

 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X2 < thr,]

 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_range <- temp_range[temp_range$X2 < thr,]

 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 temp_seasonality <- temp_seasonality[temp_seasonality$X2 < thr,]

 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))
 tmax_clim_change <- tmax_clim_change[tmax_clim_change$X2 < thr,]

 jim_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)


 res <- lapply(jim_list, function(x){
  test <- length(intersect(gab_outliers$orthogroup, x$X1))
  return(test)
 })

 # Hypergeometric expectancy
 #res_percent <- lapply(res,"/",length(annual_precip$X1))
 total <- length(unique(sort(jim$Orthogroup)))
 percent <- 100/(thr * 100)
 expect_outliers <- total/percent
 hyper_exp <- (expect_outliers * expect_outliers)/total
 res_percent <- lapply(res,"/",hyper_exp)

df11 <- data.frame(res_percent)
colnames(df11) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df11) <- c("tiffin_truncatula")






##### Preparing Gabriele top thr % orthogroup p values
###### ptremula
gab <- read.table("ingvarsson_ptremula_final_orthogroup.txt", header=TRUE)
gab_ortho_p <- gab[,7]
gab_ranks <- assign.rank(gab_ortho_p)
gab <- cbind(gab,gab_ranks)
gab_outliers <- gab[gab$gab_ranks < thr,]

##### Preparing Jim top thr % ortho p values
jim <- readRDS("run220927_COMBINED25species_updatedOF_221213_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Populus tremula",]
jim <-jim[order(jim$Orthogroup),]
jim_ranks <- sapply(split(jim$min_sdP_DS, jim$climate), assign.rank)
jim_ranks <- data.frame(jim_ranks)

 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 annual_precip <- annual_precip[annual_precip$X2 < thr,]

 isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 isothermality <- isothermality[isothermality$X2 < thr,]

 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X2 < thr,]

 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_diurnal <- mean_diurnal[mean_diurnal$X2 < thr,]

 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp <- mean_temp[mean_temp$X2 < thr,]

 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X2 < thr,]

 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X2 < thr,]

 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X2 < thr,]

 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X2 < thr,]

 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X2 < thr,]

 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 prec_clim_change <- prec_clim_change[prec_clim_change$X2 < thr,]

 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X2 < thr,]

 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_month <- precip_dry_month[precip_dry_month$X2 < thr,]

 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X2 < thr,]

 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_seasonality <- precip_seasonality[precip_seasonality$X2 < thr,]

 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X2 < thr,]

 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_month <- precip_wet_month[precip_wet_month$X2 < thr,]

 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X2 < thr,]

 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_range <- temp_range[temp_range$X2 < thr,]

 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 temp_seasonality <- temp_seasonality[temp_seasonality$X2 < thr,]

 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))
 tmax_clim_change <- tmax_clim_change[tmax_clim_change$X2 < thr,]

 jim_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)


 res <- lapply(jim_list, function(x){
  test <- length(intersect(gab_outliers$orthogroup, x$X1))
  return(test)
 })

 # Hypergeometric expectancy
 #res_percent <- lapply(res,"/",length(annual_precip$X1))
 total <- length(unique(sort(jim$Orthogroup)))
 percent <- 100/(thr * 100)
 expect_outliers <- total/percent
 hyper_exp <- (expect_outliers * expect_outliers)/total
 res_percent <- lapply(res,"/",hyper_exp)

df12 <- data.frame(res_percent)
colnames(df12) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df12) <- c("ingvarsson_ptremula")






##### Preparing Gabriele top thr % orthogroup p values
###### hannus
gab <- read.table("todesco_hann_final_orthogroup.txt", header=TRUE)
gab_ortho_p <- gab[,7]
gab_ranks <- assign.rank(gab_ortho_p)
gab <- cbind(gab,gab_ranks)
gab_outliers <- gab[gab$gab_ranks < thr,]

##### Preparing Jim top thr % ortho p values
jim <- readRDS("run220927_COMBINED25species_updatedOF_221213_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Helianthus annuus",]
jim <-jim[order(jim$Orthogroup),]
jim_ranks <- sapply(split(jim$min_sdP_DS, jim$climate), assign.rank)
jim_ranks <- data.frame(jim_ranks)

 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 annual_precip <- annual_precip[annual_precip$X2 < thr,]

 isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 isothermality <- isothermality[isothermality$X2 < thr,]

 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X2 < thr,]

 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_diurnal <- mean_diurnal[mean_diurnal$X2 < thr,]

 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp <- mean_temp[mean_temp$X2 < thr,]

 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X2 < thr,]

 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X2 < thr,]

 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X2 < thr,]

 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X2 < thr,]

 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X2 < thr,]

 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 prec_clim_change <- prec_clim_change[prec_clim_change$X2 < thr,]

 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X2 < thr,]

 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_month <- precip_dry_month[precip_dry_month$X2 < thr,]

 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X2 < thr,]

 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_seasonality <- precip_seasonality[precip_seasonality$X2 < thr,]

 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X2 < thr,]

 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_month <- precip_wet_month[precip_wet_month$X2 < thr,]

 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X2 < thr,]

 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_range <- temp_range[temp_range$X2 < thr,]

 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 temp_seasonality <- temp_seasonality[temp_seasonality$X2 < thr,]

 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))
 tmax_clim_change <- tmax_clim_change[tmax_clim_change$X2 < thr,]

 jim_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)


 res <- lapply(jim_list, function(x){
  test <- length(intersect(gab_outliers$orthogroup, x$X1))
  return(test)
 })

 # Hypergeometric expectancy
 #res_percent <- lapply(res,"/",length(annual_precip$X1))
 total <- length(unique(sort(jim$Orthogroup)))
 percent <- 100/(thr * 100)
 expect_outliers <- total/percent
 hyper_exp <- (expect_outliers * expect_outliers)/total
 res_percent <- lapply(res,"/",hyper_exp)

df13 <- data.frame(res_percent)
colnames(df13) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df13) <- c("todesco_hann")







##### Preparing Gabriele top thr % orthogroup p values
###### hannus
gab <- read.table("todesco_hpet_final_orthogroup.txt", header=TRUE)
gab_ortho_p <- gab[,7]
gab_ranks <- assign.rank(gab_ortho_p)
gab <- cbind(gab,gab_ranks)
gab_outliers <- gab[gab$gab_ranks < thr,]

##### Preparing Jim top thr % ortho p values
jim <- readRDS("run220927_COMBINED25species_updatedOF_221213_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Helianthus petiolaris",]
jim <-jim[order(jim$Orthogroup),]
jim_ranks <- sapply(split(jim$min_sdP_DS, jim$climate), assign.rank)
jim_ranks <- data.frame(jim_ranks)

 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 annual_precip <- annual_precip[annual_precip$X2 < thr,]

 isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 isothermality <- isothermality[isothermality$X2 < thr,]

 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X2 < thr,]

 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_diurnal <- mean_diurnal[mean_diurnal$X2 < thr,]

 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp <- mean_temp[mean_temp$X2 < thr,]

 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X2 < thr,]

 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X2 < thr,]

 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X2 < thr,]

 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X2 < thr,]

 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X2 < thr,]

 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 prec_clim_change <- prec_clim_change[prec_clim_change$X2 < thr,]

 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X2 < thr,]

 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_month <- precip_dry_month[precip_dry_month$X2 < thr,]

 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X2 < thr,]

 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_seasonality <- precip_seasonality[precip_seasonality$X2 < thr,]

 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X2 < thr,]

 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_month <- precip_wet_month[precip_wet_month$X2 < thr,]

 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X2 < thr,]

 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_range <- temp_range[temp_range$X2 < thr,]

 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 temp_seasonality <- temp_seasonality[temp_seasonality$X2 < thr,]

 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))
 tmax_clim_change <- tmax_clim_change[tmax_clim_change$X2 < thr,]

 jim_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)


 res <- lapply(jim_list, function(x){
  test <- length(intersect(gab_outliers$orthogroup, x$X1))
  return(test)
 })

 # Hypergeometric expectancy
 #res_percent <- lapply(res,"/",length(annual_precip$X1))
 total <- length(unique(sort(jim$Orthogroup)))
 percent <- 100/(thr * 100)
 expect_outliers <- total/percent
 hyper_exp <- (expect_outliers * expect_outliers)/total
 res_percent <- lapply(res,"/",hyper_exp)

df14 <- data.frame(res_percent)
colnames(df14) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df14) <- c("todesco_hpet")






##### Preparing Gabriele top thr % orthogroup p values
###### hannus
gab <- read.table("todesco_hargo_final_orthogroup.txt", header=TRUE)
gab_ortho_p <- gab[,7]
gab_ranks <- assign.rank(gab_ortho_p)
gab <- cbind(gab,gab_ranks)
gab_outliers <- gab[gab$gab_ranks < thr,]

##### Preparing Jim top thr % ortho p values
jim <- readRDS("run220927_COMBINED25species_updatedOF_221213_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Helianthus argophyllus",]
jim <-jim[order(jim$Orthogroup),]
jim_ranks <- sapply(split(jim$min_sdP_DS, jim$climate), assign.rank)
jim_ranks <- data.frame(jim_ranks)

 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 annual_precip <- annual_precip[annual_precip$X2 < thr,]

 isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 isothermality <- isothermality[isothermality$X2 < thr,]

 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X2 < thr,]

 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_diurnal <- mean_diurnal[mean_diurnal$X2 < thr,]

 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp <- mean_temp[mean_temp$X2 < thr,]

 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X2 < thr,]

 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X2 < thr,]

 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X2 < thr,]

 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X2 < thr,]

 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X2 < thr,]

 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 prec_clim_change <- prec_clim_change[prec_clim_change$X2 < thr,]

 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X2 < thr,]

 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_month <- precip_dry_month[precip_dry_month$X2 < thr,]

 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X2 < thr,]

 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_seasonality <- precip_seasonality[precip_seasonality$X2 < thr,]

 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X2 < thr,]

 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_month <- precip_wet_month[precip_wet_month$X2 < thr,]

 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X2 < thr,]

 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_range <- temp_range[temp_range$X2 < thr,]

 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 temp_seasonality <- temp_seasonality[temp_seasonality$X2 < thr,]

 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))
 tmax_clim_change <- tmax_clim_change[tmax_clim_change$X2 < thr,]

 jim_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)


 res <- lapply(jim_list, function(x){
  test <- length(intersect(gab_outliers$orthogroup, x$X1))
  return(test)
 })

 # Hypergeometric expectancy
 #res_percent <- lapply(res,"/",length(annual_precip$X1))
 total <- length(unique(sort(jim$Orthogroup)))
 percent <- 100/(thr * 100)
 expect_outliers <- total/percent
 hyper_exp <- (expect_outliers * expect_outliers)/total
 res_percent <- lapply(res,"/",hyper_exp)

df15 <- data.frame(res_percent)
colnames(df15) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df15) <- c("todesco_hargo")





# df8 removed, no results for alyrata from SweeD
final <- rbind(df,df2,df3,df4,df5,df6,df7,df9,df10,df11,df12,df13,df14,df15)
final <- as.matrix(final)


Heatmap(final, show_column_dend = FALSE,show_row_dend = FALSE,heatmap_legend_param = list(title = "", color_bar = "continuous"),row_names_side = "left")
