

### Preparing Gabriele top 10% p values

library(ComplexHeatmap)
gab <- read.table("final_orthogroup.txt", header=TRUE)

gab_ortho_p <- gab[,7]
thr <- 0.1

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


gab_ranks <- assign.rank(gab_ortho_p)

gab <- cbind(gab,gab_ranks)


gab_outliers <- gab[gab$gab_ranks < thr,]


### Preparing Jim top 10% p values
jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
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

 res_percent <- lapply(res,"/",length(annual_precip$X1))

df <- data.frame(res_percent)
colnames(df) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df) <- c("wright_atuber")


#########################################
jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
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

 res_percent <- lapply(res,"/",length(annual_precip$X1))

df2 <- data.frame(res_percent)
colnames(df2) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df2) <- c("lowry_phallii")


 #####################################
 jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
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

  res_percent <- lapply(res,"/",length(annual_precip$X1))

 df3 <- data.frame(res_percent)
 colnames(df3) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
 rownames(df3) <- c("weigel_athaliana")

########

jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
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

 res_percent <- lapply(res,"/",length(annual_precip$X1))

df4 <- data.frame(res_percent)
colnames(df4) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df4) <- c("weigel_capsella")


######
jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
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

 res_percent <- lapply(res,"/",length(annual_precip$X1))

df5 <- data.frame(res_percent)
colnames(df5) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df5) <- c("murray_emag")
#########


jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
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

 res_percent <- lapply(res,"/",length(annual_precip$X1))

df6 <- data.frame(res_percent)
colnames(df6) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df6) <- c("kubota_ahalleri")


######
jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
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

 res_percent <- lapply(res,"/",length(annual_precip$X1))

df7 <- data.frame(res_percent)
colnames(df7) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df7) <- c("murray_esid")



#####


jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Arabidopsis lyrata",]

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

 res_percent <- lapply(res,"/",length(annual_precip$X1))

df8 <- data.frame(res_percent)
colnames(df8) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df8) <- c("savolainen_alyrata")


####

jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
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

 res_percent <- lapply(res,"/",length(annual_precip$X1))

df9 <- data.frame(res_percent)
colnames(df9) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df9) <- c("mitchell_bstricta")


#####

jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
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

 res_percent <- lapply(res,"/",length(annual_precip$X1))

df10 <- data.frame(res_percent)
colnames(df10) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df10) <- c("murray_ealb")


#####
jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
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

 res_percent <- lapply(res,"/",length(annual_precip$X1))

df11 <- data.frame(res_percent)
colnames(df11) <- c("annual_precip","isothermality", "max_temp_warmest_month", "mean_diurnal", "mean_temp", "mean_temp_cold_quarter", "mean_temp_dry_quarter", "mean_temp_warm_quarter", "mean_temp_wet_quarter", "min_temp_coldest_month", "prec_clim_change","precip_cold_quarter","precip_dry_month","precip_dry_quarter","precip_seasonality","precip_warm_quarter","precip_wet_month","precip_wet_quarter","temp_range","temp_seasonality","tmax_clim_change")
rownames(df11) <- c("tiffin_truncatula")



final <- rbind(df,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11)
final <- as.matrix(final)

Heatmap(final, show_column_dend = FALSE,show_row_dend = FALSE,heatmap_legend_param = list(title = "", color_bar = "discrete"))
