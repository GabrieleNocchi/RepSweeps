    ##### change wd accordingly
    library(dplyr)
    setwd("/Users/gnocc/Desktop/projects/jim_results_and_ortho_map/")

    ##### Preparing Gabriele orthogroup p values
    ##### wright_atuber
    gab <- read.table("wright_atuber_final_orthogroup.txt", header=TRUE)




    ##### Preparing Jim ortho p values
    jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
    jim <- jim[jim$species == "Amaranthus tuberculatus",]

    jim <-jim[order(jim$Orthogroup),]

    jim_ranks <- split(jim$min_sdP_DS, jim$climate)

    jim_ranks <- data.frame(jim_ranks)




     annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
     isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
     max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
     mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
     mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
     mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
     mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
     mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
     mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
     min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
     prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
     precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
     precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
     precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
     precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
     precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
     precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
     precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
     temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
     temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
     tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))




      test <- intersect(gab$orthogroup, annual_precip$X1)
      annual_precip <- annual_precip[annual_precip$X1 %in% test,]
      colnames(annual_precip) <- c("orthogroup", "min_sdP_DS")
      annual_precip <- annual_precip %>%
      left_join(gab, by = c("orthogroup"))

      test <- intersect(gab$orthogroup, isothermality$X1)
      isothermality <- isothermality[isothermality$X1 %in% test,]
      colnames(isothermality) <- c("orthogroup", "min_sdP_DS")
      isothermality <- isothermality %>%
      left_join(gab, by = c("orthogroup"))

      test <- intersect(gab$orthogroup, max_temp_warmest_month$X1)
      max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X1 %in% test,]
      colnames(max_temp_warmest_month) <- c("orthogroup", "min_sdP_DS")
      max_temp_warmest_month <- max_temp_warmest_month %>%
      left_join(gab, by = c("orthogroup"))


      test <- intersect(gab$orthogroup, mean_diurnal$X1)
      mean_diurnal <- mean_diurnal[mean_diurnal$X1 %in% test,]
      colnames(mean_diurnal) <- c("orthogroup", "min_sdP_DS")
      mean_diurnal <- mean_diurnal %>%
      left_join(gab, by = c("orthogroup"))


      test <- intersect(gab$orthogroup, mean_temp$X1)
      mean_temp <- mean_temp[mean_temp$X1 %in% test,]
      colnames(mean_temp) <- c("orthogroup", "min_sdP_DS")
      mean_temp <- mean_temp %>%
      left_join(gab, by = c("orthogroup"))


      test <- intersect(gab$orthogroup, mean_temp_cold_quarter$X1)
      mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X1 %in% test,]
      colnames(mean_temp_cold_quarter) <- c("orthogroup", "min_sdP_DS")
      mean_temp_cold_quarter <- mean_temp_cold_quarter %>%
      left_join(gab, by = c("orthogroup"))


      test <- intersect(gab$orthogroup, mean_temp_dry_quarter$X1)
      mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X1 %in% test,]
      colnames(mean_temp_dry_quarter) <- c("orthogroup", "min_sdP_DS")
      mean_temp_dry_quarter <- mean_temp_dry_quarter %>%
      left_join(gab, by = c("orthogroup"))


      test <- intersect(gab$orthogroup, mean_temp_warm_quarter$X1)
      mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X1 %in% test,]
      colnames(mean_temp_warm_quarter) <- c("orthogroup", "min_sdP_DS")
      mean_temp_warm_quarter <- mean_temp_warm_quarter %>%
      left_join(gab, by = c("orthogroup"))


      test <- intersect(gab$orthogroup, mean_temp_wet_quarter$X1)
      mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X1 %in% test,]
      colnames(mean_temp_wet_quarter) <- c("orthogroup", "min_sdP_DS")
      mean_temp_wet_quarter <- mean_temp_wet_quarter %>%
      left_join(gab, by = c("orthogroup"))


      test <- intersect(gab$orthogroup, min_temp_coldest_month$X1)
      min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X1 %in% test,]
      colnames(min_temp_coldest_month) <- c("orthogroup", "min_sdP_DS")
      min_temp_coldest_month <- min_temp_coldest_month %>%
      left_join(gab, by = c("orthogroup"))


      test <- intersect(gab$orthogroup, prec_clim_change$X1)
      prec_clim_change <- prec_clim_change[prec_clim_change$X1 %in% test,]
      colnames(prec_clim_change) <- c("orthogroup", "min_sdP_DS")
      prec_clim_change <- prec_clim_change %>%
      left_join(gab, by = c("orthogroup"))


      test <- intersect(gab$orthogroup, precip_cold_quarter$X1)
      precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X1 %in% test,]
      colnames(precip_cold_quarter) <- c("orthogroup", "min_sdP_DS")
      precip_cold_quarter <- precip_cold_quarter %>%
      left_join(gab, by = c("orthogroup"))


      test <- intersect(gab$orthogroup, precip_dry_month$X1)
      precip_dry_month <- precip_dry_month[precip_dry_month$X1 %in% test,]
      colnames(precip_dry_month) <- c("orthogroup", "min_sdP_DS")
      precip_dry_month <- precip_dry_month %>%
      left_join(gab, by = c("orthogroup"))


      test <- intersect(gab$orthogroup, precip_dry_quarter$X1)
      precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X1 %in% test,]
      colnames(precip_dry_quarter) <- c("orthogroup", "min_sdP_DS")
      precip_dry_quarter <- precip_dry_quarter %>%
      left_join(gab, by = c("orthogroup"))


      test <- intersect(gab$orthogroup, precip_seasonality$X1)
      precip_seasonality <- precip_seasonality[precip_seasonality$X1 %in% test,]
      colnames(precip_seasonality) <- c("orthogroup", "min_sdP_DS")
      precip_seasonality <- precip_seasonality %>%
      left_join(gab, by = c("orthogroup"))



      test <- intersect(gab$orthogroup, precip_warm_quarter$X1)
      precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X1 %in% test,]
      colnames(precip_warm_quarter) <- c("orthogroup", "min_sdP_DS")
      precip_warm_quarter <- precip_warm_quarter %>%
      left_join(gab, by = c("orthogroup"))



      test <- intersect(gab$orthogroup, precip_wet_month$X1)
      precip_wet_month <- precip_wet_month[precip_wet_month$X1 %in% test,]
      colnames(precip_wet_month) <- c("orthogroup", "min_sdP_DS")
      precip_wet_month <- precip_wet_month %>%
      left_join(gab, by = c("orthogroup"))


      test <- intersect(gab$orthogroup, precip_wet_quarter$X1)
      precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X1 %in% test,]
      colnames(precip_wet_quarter) <- c("orthogroup", "min_sdP_DS")
      precip_wet_quarter <- precip_wet_quarter %>%
      left_join(gab, by = c("orthogroup"))


      test <- intersect(gab$orthogroup, temp_range$X1)
      temp_range <- temp_range[temp_range$X1 %in% test,]
      colnames(temp_range) <- c("orthogroup", "min_sdP_DS")
      temp_range <- temp_range %>%
      left_join(gab, by = c("orthogroup"))



      test <- intersect(gab$orthogroup, temp_seasonality$X1)
      temp_seasonality <- temp_seasonality[temp_seasonality$X1 %in% test,]
      colnames(temp_seasonality) <- c("orthogroup", "min_sdP_DS")
      temp_seasonality <- temp_seasonality %>%
      left_join(gab, by = c("orthogroup"))


      test <- intersect(gab$orthogroup, tmax_clim_change$X1)
      tmax_clim_change <- tmax_clim_change[tmax_clim_change$X1 %in% test,]
      colnames(tmax_clim_change) <- c("orthogroup", "min_sdP_DS")
      tmax_clim_change <- tmax_clim_change %>%
      left_join(gab, by = c("orthogroup"))


    my_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)

    pdf("wright_atuber_scatter.pdf", height = 14, width = 14)
  par(mfrow = c(7, 3), mar=c(4,4,4,4), cex.axis=0.7, cex.lab=1, cex.main=1, cex.sub=1,mgp=c(2,1,0) )
  lapply(my_list, function(x){
      return(plot(-log10(as.numeric(x$min_sdP_DS)),-log10(x$dunnsidak_orthopvalues), xlab="Jim_WZA_ortho_p", ylab="Gabriele_SweeD_ortho_p"))
  })

  dev.off()



  ##### Preparing Gabriele orthogroup p values
  ##### lowry_phalli
  gab <- read.table("lowry_phalli_final_orthogroup.txt", header=TRUE)




  ##### Preparing Jim ortho p values
  jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
  jim <- jim[jim$species == "Panicum hallii",]

  jim <-jim[order(jim$Orthogroup),]

  jim_ranks <- split(jim$min_sdP_DS, jim$climate)

  jim_ranks <- data.frame(jim_ranks)




   annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
   isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
   max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
   mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
   mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
   mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
   mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
   mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
   mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
   min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
   prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
   precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
   precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
   precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
   precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
   precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
   precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
   precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
   temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
   temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
   tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))




    test <- intersect(gab$orthogroup, annual_precip$X1)
    annual_precip <- annual_precip[annual_precip$X1 %in% test,]
    colnames(annual_precip) <- c("orthogroup", "min_sdP_DS")
    annual_precip <- annual_precip %>%
    left_join(gab, by = c("orthogroup"))

    test <- intersect(gab$orthogroup, isothermality$X1)
    isothermality <- isothermality[isothermality$X1 %in% test,]
    colnames(isothermality) <- c("orthogroup", "min_sdP_DS")
    isothermality <- isothermality %>%
    left_join(gab, by = c("orthogroup"))

    test <- intersect(gab$orthogroup, max_temp_warmest_month$X1)
    max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X1 %in% test,]
    colnames(max_temp_warmest_month) <- c("orthogroup", "min_sdP_DS")
    max_temp_warmest_month <- max_temp_warmest_month %>%
    left_join(gab, by = c("orthogroup"))


    test <- intersect(gab$orthogroup, mean_diurnal$X1)
    mean_diurnal <- mean_diurnal[mean_diurnal$X1 %in% test,]
    colnames(mean_diurnal) <- c("orthogroup", "min_sdP_DS")
    mean_diurnal <- mean_diurnal %>%
    left_join(gab, by = c("orthogroup"))


    test <- intersect(gab$orthogroup, mean_temp$X1)
    mean_temp <- mean_temp[mean_temp$X1 %in% test,]
    colnames(mean_temp) <- c("orthogroup", "min_sdP_DS")
    mean_temp <- mean_temp %>%
    left_join(gab, by = c("orthogroup"))


    test <- intersect(gab$orthogroup, mean_temp_cold_quarter$X1)
    mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X1 %in% test,]
    colnames(mean_temp_cold_quarter) <- c("orthogroup", "min_sdP_DS")
    mean_temp_cold_quarter <- mean_temp_cold_quarter %>%
    left_join(gab, by = c("orthogroup"))


    test <- intersect(gab$orthogroup, mean_temp_dry_quarter$X1)
    mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X1 %in% test,]
    colnames(mean_temp_dry_quarter) <- c("orthogroup", "min_sdP_DS")
    mean_temp_dry_quarter <- mean_temp_dry_quarter %>%
    left_join(gab, by = c("orthogroup"))


    test <- intersect(gab$orthogroup, mean_temp_warm_quarter$X1)
    mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X1 %in% test,]
    colnames(mean_temp_warm_quarter) <- c("orthogroup", "min_sdP_DS")
    mean_temp_warm_quarter <- mean_temp_warm_quarter %>%
    left_join(gab, by = c("orthogroup"))


    test <- intersect(gab$orthogroup, mean_temp_wet_quarter$X1)
    mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X1 %in% test,]
    colnames(mean_temp_wet_quarter) <- c("orthogroup", "min_sdP_DS")
    mean_temp_wet_quarter <- mean_temp_wet_quarter %>%
    left_join(gab, by = c("orthogroup"))


    test <- intersect(gab$orthogroup, min_temp_coldest_month$X1)
    min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X1 %in% test,]
    colnames(min_temp_coldest_month) <- c("orthogroup", "min_sdP_DS")
    min_temp_coldest_month <- min_temp_coldest_month %>%
    left_join(gab, by = c("orthogroup"))


    test <- intersect(gab$orthogroup, prec_clim_change$X1)
    prec_clim_change <- prec_clim_change[prec_clim_change$X1 %in% test,]
    colnames(prec_clim_change) <- c("orthogroup", "min_sdP_DS")
    prec_clim_change <- prec_clim_change %>%
    left_join(gab, by = c("orthogroup"))


    test <- intersect(gab$orthogroup, precip_cold_quarter$X1)
    precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X1 %in% test,]
    colnames(precip_cold_quarter) <- c("orthogroup", "min_sdP_DS")
    precip_cold_quarter <- precip_cold_quarter %>%
    left_join(gab, by = c("orthogroup"))


    test <- intersect(gab$orthogroup, precip_dry_month$X1)
    precip_dry_month <- precip_dry_month[precip_dry_month$X1 %in% test,]
    colnames(precip_dry_month) <- c("orthogroup", "min_sdP_DS")
    precip_dry_month <- precip_dry_month %>%
    left_join(gab, by = c("orthogroup"))


    test <- intersect(gab$orthogroup, precip_dry_quarter$X1)
    precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X1 %in% test,]
    colnames(precip_dry_quarter) <- c("orthogroup", "min_sdP_DS")
    precip_dry_quarter <- precip_dry_quarter %>%
    left_join(gab, by = c("orthogroup"))


    test <- intersect(gab$orthogroup, precip_seasonality$X1)
    precip_seasonality <- precip_seasonality[precip_seasonality$X1 %in% test,]
    colnames(precip_seasonality) <- c("orthogroup", "min_sdP_DS")
    precip_seasonality <- precip_seasonality %>%
    left_join(gab, by = c("orthogroup"))



    test <- intersect(gab$orthogroup, precip_warm_quarter$X1)
    precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X1 %in% test,]
    colnames(precip_warm_quarter) <- c("orthogroup", "min_sdP_DS")
    precip_warm_quarter <- precip_warm_quarter %>%
    left_join(gab, by = c("orthogroup"))



    test <- intersect(gab$orthogroup, precip_wet_month$X1)
    precip_wet_month <- precip_wet_month[precip_wet_month$X1 %in% test,]
    colnames(precip_wet_month) <- c("orthogroup", "min_sdP_DS")
    precip_wet_month <- precip_wet_month %>%
    left_join(gab, by = c("orthogroup"))


    test <- intersect(gab$orthogroup, precip_wet_quarter$X1)
    precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X1 %in% test,]
    colnames(precip_wet_quarter) <- c("orthogroup", "min_sdP_DS")
    precip_wet_quarter <- precip_wet_quarter %>%
    left_join(gab, by = c("orthogroup"))


    test <- intersect(gab$orthogroup, temp_range$X1)
    temp_range <- temp_range[temp_range$X1 %in% test,]
    colnames(temp_range) <- c("orthogroup", "min_sdP_DS")
    temp_range <- temp_range %>%
    left_join(gab, by = c("orthogroup"))



    test <- intersect(gab$orthogroup, temp_seasonality$X1)
    temp_seasonality <- temp_seasonality[temp_seasonality$X1 %in% test,]
    colnames(temp_seasonality) <- c("orthogroup", "min_sdP_DS")
    temp_seasonality <- temp_seasonality %>%
    left_join(gab, by = c("orthogroup"))


    test <- intersect(gab$orthogroup, tmax_clim_change$X1)
    tmax_clim_change <- tmax_clim_change[tmax_clim_change$X1 %in% test,]
    colnames(tmax_clim_change) <- c("orthogroup", "min_sdP_DS")
    tmax_clim_change <- tmax_clim_change %>%
    left_join(gab, by = c("orthogroup"))


  my_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)

  pdf("lowry_phalli.pdf", height = 14, width = 14)
  par(mfrow = c(7, 3), mar=c(4,4,4,4), cex.axis=0.7, cex.lab=1, cex.main=1, cex.sub=1,mgp=c(2,1,0) )
  lapply(my_list, function(x){
      return(plot(-log10(as.numeric(x$min_sdP_DS)),-log10(x$dunnsidak_orthopvalues), xlab="Jim_WZA_ortho_p", ylab="Gabriele_SweeD_ortho_p"))
  })

  dev.off()



##### Preparing Gabriele orthogroup p values
##### weigel_athal
gab <- read.table("weigel_athaliana_IBE_final_orthogroup.txt", header=TRUE)




##### Preparing Jim ortho p values
jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Arabidopsis thaliana",]

jim <-jim[order(jim$Orthogroup),]

jim_ranks <- split(jim$min_sdP_DS, jim$climate)

jim_ranks <- data.frame(jim_ranks)




 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))




  test <- intersect(gab$orthogroup, annual_precip$X1)
  annual_precip <- annual_precip[annual_precip$X1 %in% test,]
  colnames(annual_precip) <- c("orthogroup", "min_sdP_DS")
  annual_precip <- annual_precip %>%
  left_join(gab, by = c("orthogroup"))

  test <- intersect(gab$orthogroup, isothermality$X1)
  isothermality <- isothermality[isothermality$X1 %in% test,]
  colnames(isothermality) <- c("orthogroup", "min_sdP_DS")
  isothermality <- isothermality %>%
  left_join(gab, by = c("orthogroup"))

  test <- intersect(gab$orthogroup, max_temp_warmest_month$X1)
  max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X1 %in% test,]
  colnames(max_temp_warmest_month) <- c("orthogroup", "min_sdP_DS")
  max_temp_warmest_month <- max_temp_warmest_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_diurnal$X1)
  mean_diurnal <- mean_diurnal[mean_diurnal$X1 %in% test,]
  colnames(mean_diurnal) <- c("orthogroup", "min_sdP_DS")
  mean_diurnal <- mean_diurnal %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp$X1)
  mean_temp <- mean_temp[mean_temp$X1 %in% test,]
  colnames(mean_temp) <- c("orthogroup", "min_sdP_DS")
  mean_temp <- mean_temp %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_cold_quarter$X1)
  mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X1 %in% test,]
  colnames(mean_temp_cold_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_cold_quarter <- mean_temp_cold_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_dry_quarter$X1)
  mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X1 %in% test,]
  colnames(mean_temp_dry_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_dry_quarter <- mean_temp_dry_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_warm_quarter$X1)
  mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X1 %in% test,]
  colnames(mean_temp_warm_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_warm_quarter <- mean_temp_warm_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_wet_quarter$X1)
  mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X1 %in% test,]
  colnames(mean_temp_wet_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_wet_quarter <- mean_temp_wet_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, min_temp_coldest_month$X1)
  min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X1 %in% test,]
  colnames(min_temp_coldest_month) <- c("orthogroup", "min_sdP_DS")
  min_temp_coldest_month <- min_temp_coldest_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, prec_clim_change$X1)
  prec_clim_change <- prec_clim_change[prec_clim_change$X1 %in% test,]
  colnames(prec_clim_change) <- c("orthogroup", "min_sdP_DS")
  prec_clim_change <- prec_clim_change %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_cold_quarter$X1)
  precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X1 %in% test,]
  colnames(precip_cold_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_cold_quarter <- precip_cold_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_dry_month$X1)
  precip_dry_month <- precip_dry_month[precip_dry_month$X1 %in% test,]
  colnames(precip_dry_month) <- c("orthogroup", "min_sdP_DS")
  precip_dry_month <- precip_dry_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_dry_quarter$X1)
  precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X1 %in% test,]
  colnames(precip_dry_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_dry_quarter <- precip_dry_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_seasonality$X1)
  precip_seasonality <- precip_seasonality[precip_seasonality$X1 %in% test,]
  colnames(precip_seasonality) <- c("orthogroup", "min_sdP_DS")
  precip_seasonality <- precip_seasonality %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, precip_warm_quarter$X1)
  precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X1 %in% test,]
  colnames(precip_warm_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_warm_quarter <- precip_warm_quarter %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, precip_wet_month$X1)
  precip_wet_month <- precip_wet_month[precip_wet_month$X1 %in% test,]
  colnames(precip_wet_month) <- c("orthogroup", "min_sdP_DS")
  precip_wet_month <- precip_wet_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_wet_quarter$X1)
  precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X1 %in% test,]
  colnames(precip_wet_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_wet_quarter <- precip_wet_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, temp_range$X1)
  temp_range <- temp_range[temp_range$X1 %in% test,]
  colnames(temp_range) <- c("orthogroup", "min_sdP_DS")
  temp_range <- temp_range %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, temp_seasonality$X1)
  temp_seasonality <- temp_seasonality[temp_seasonality$X1 %in% test,]
  colnames(temp_seasonality) <- c("orthogroup", "min_sdP_DS")
  temp_seasonality <- temp_seasonality %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, tmax_clim_change$X1)
  tmax_clim_change <- tmax_clim_change[tmax_clim_change$X1 %in% test,]
  colnames(tmax_clim_change) <- c("orthogroup", "min_sdP_DS")
  tmax_clim_change <- tmax_clim_change %>%
  left_join(gab, by = c("orthogroup"))


my_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)

pdf("weigel_athaliana_IBE_scatter.pdf", height = 14, width = 14)
par(mfrow = c(7, 3), mar=c(4,4,4,4), cex.axis=0.7, cex.lab=1, cex.main=1, cex.sub=1,mgp=c(2,1,0) )
lapply(my_list, function(x){
    return(plot(-log10(as.numeric(x$min_sdP_DS)),-log10(x$dunnsidak_orthopvalues), xlab="Jim_WZA_ortho_p", ylab="Gabriele_SweeD_ortho_p"))
})

dev.off()


##### Preparing Gabriele orthogroup p values
##### weiger_capsella
gab <- read.table("weiger_capsella_final_orthogroup.txt", header=TRUE)




##### Preparing Jim ortho p values
jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Capsella rubella",]

jim <-jim[order(jim$Orthogroup),]

jim_ranks <- split(jim$min_sdP_DS, jim$climate)

jim_ranks <- data.frame(jim_ranks)




 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))




  test <- intersect(gab$orthogroup, annual_precip$X1)
  annual_precip <- annual_precip[annual_precip$X1 %in% test,]
  colnames(annual_precip) <- c("orthogroup", "min_sdP_DS")
  annual_precip <- annual_precip %>%
  left_join(gab, by = c("orthogroup"))

  test <- intersect(gab$orthogroup, isothermality$X1)
  isothermality <- isothermality[isothermality$X1 %in% test,]
  colnames(isothermality) <- c("orthogroup", "min_sdP_DS")
  isothermality <- isothermality %>%
  left_join(gab, by = c("orthogroup"))

  test <- intersect(gab$orthogroup, max_temp_warmest_month$X1)
  max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X1 %in% test,]
  colnames(max_temp_warmest_month) <- c("orthogroup", "min_sdP_DS")
  max_temp_warmest_month <- max_temp_warmest_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_diurnal$X1)
  mean_diurnal <- mean_diurnal[mean_diurnal$X1 %in% test,]
  colnames(mean_diurnal) <- c("orthogroup", "min_sdP_DS")
  mean_diurnal <- mean_diurnal %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp$X1)
  mean_temp <- mean_temp[mean_temp$X1 %in% test,]
  colnames(mean_temp) <- c("orthogroup", "min_sdP_DS")
  mean_temp <- mean_temp %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_cold_quarter$X1)
  mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X1 %in% test,]
  colnames(mean_temp_cold_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_cold_quarter <- mean_temp_cold_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_dry_quarter$X1)
  mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X1 %in% test,]
  colnames(mean_temp_dry_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_dry_quarter <- mean_temp_dry_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_warm_quarter$X1)
  mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X1 %in% test,]
  colnames(mean_temp_warm_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_warm_quarter <- mean_temp_warm_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_wet_quarter$X1)
  mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X1 %in% test,]
  colnames(mean_temp_wet_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_wet_quarter <- mean_temp_wet_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, min_temp_coldest_month$X1)
  min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X1 %in% test,]
  colnames(min_temp_coldest_month) <- c("orthogroup", "min_sdP_DS")
  min_temp_coldest_month <- min_temp_coldest_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, prec_clim_change$X1)
  prec_clim_change <- prec_clim_change[prec_clim_change$X1 %in% test,]
  colnames(prec_clim_change) <- c("orthogroup", "min_sdP_DS")
  prec_clim_change <- prec_clim_change %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_cold_quarter$X1)
  precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X1 %in% test,]
  colnames(precip_cold_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_cold_quarter <- precip_cold_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_dry_month$X1)
  precip_dry_month <- precip_dry_month[precip_dry_month$X1 %in% test,]
  colnames(precip_dry_month) <- c("orthogroup", "min_sdP_DS")
  precip_dry_month <- precip_dry_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_dry_quarter$X1)
  precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X1 %in% test,]
  colnames(precip_dry_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_dry_quarter <- precip_dry_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_seasonality$X1)
  precip_seasonality <- precip_seasonality[precip_seasonality$X1 %in% test,]
  colnames(precip_seasonality) <- c("orthogroup", "min_sdP_DS")
  precip_seasonality <- precip_seasonality %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, precip_warm_quarter$X1)
  precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X1 %in% test,]
  colnames(precip_warm_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_warm_quarter <- precip_warm_quarter %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, precip_wet_month$X1)
  precip_wet_month <- precip_wet_month[precip_wet_month$X1 %in% test,]
  colnames(precip_wet_month) <- c("orthogroup", "min_sdP_DS")
  precip_wet_month <- precip_wet_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_wet_quarter$X1)
  precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X1 %in% test,]
  colnames(precip_wet_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_wet_quarter <- precip_wet_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, temp_range$X1)
  temp_range <- temp_range[temp_range$X1 %in% test,]
  colnames(temp_range) <- c("orthogroup", "min_sdP_DS")
  temp_range <- temp_range %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, temp_seasonality$X1)
  temp_seasonality <- temp_seasonality[temp_seasonality$X1 %in% test,]
  colnames(temp_seasonality) <- c("orthogroup", "min_sdP_DS")
  temp_seasonality <- temp_seasonality %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, tmax_clim_change$X1)
  tmax_clim_change <- tmax_clim_change[tmax_clim_change$X1 %in% test,]
  colnames(tmax_clim_change) <- c("orthogroup", "min_sdP_DS")
  tmax_clim_change <- tmax_clim_change %>%
  left_join(gab, by = c("orthogroup"))


my_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)

pdf("capsella_rubella_scatter.pdf", height = 14, width = 14)
par(mfrow = c(7, 3), mar=c(4,4,4,4), cex.axis=0.7, cex.lab=1, cex.main=1, cex.sub=1,mgp=c(2,1,0) )
lapply(my_list, function(x){
    return(plot(-log10(as.numeric(x$min_sdP_DS)),-log10(x$dunnsidak_orthopvalues), xlab="Jim_WZA_ortho_p", ylab="Gabriele_SweeD_ortho_p"))
})

dev.off()


##### Preparing Gabriele orthogroup p values
##### kubota_ahalleri
gab <- read.table("kubota_ahalleri_final_orthogroup.txt", header=TRUE)




##### Preparing Jim ortho p values
jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Arabidopsis halleri",]

jim <-jim[order(jim$Orthogroup),]

jim_ranks <- split(jim$min_sdP_DS, jim$climate)

jim_ranks <- data.frame(jim_ranks)




 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))




  test <- intersect(gab$orthogroup, annual_precip$X1)
  annual_precip <- annual_precip[annual_precip$X1 %in% test,]
  colnames(annual_precip) <- c("orthogroup", "min_sdP_DS")
  annual_precip <- annual_precip %>%
  left_join(gab, by = c("orthogroup"))

  test <- intersect(gab$orthogroup, isothermality$X1)
  isothermality <- isothermality[isothermality$X1 %in% test,]
  colnames(isothermality) <- c("orthogroup", "min_sdP_DS")
  isothermality <- isothermality %>%
  left_join(gab, by = c("orthogroup"))

  test <- intersect(gab$orthogroup, max_temp_warmest_month$X1)
  max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X1 %in% test,]
  colnames(max_temp_warmest_month) <- c("orthogroup", "min_sdP_DS")
  max_temp_warmest_month <- max_temp_warmest_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_diurnal$X1)
  mean_diurnal <- mean_diurnal[mean_diurnal$X1 %in% test,]
  colnames(mean_diurnal) <- c("orthogroup", "min_sdP_DS")
  mean_diurnal <- mean_diurnal %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp$X1)
  mean_temp <- mean_temp[mean_temp$X1 %in% test,]
  colnames(mean_temp) <- c("orthogroup", "min_sdP_DS")
  mean_temp <- mean_temp %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_cold_quarter$X1)
  mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X1 %in% test,]
  colnames(mean_temp_cold_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_cold_quarter <- mean_temp_cold_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_dry_quarter$X1)
  mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X1 %in% test,]
  colnames(mean_temp_dry_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_dry_quarter <- mean_temp_dry_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_warm_quarter$X1)
  mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X1 %in% test,]
  colnames(mean_temp_warm_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_warm_quarter <- mean_temp_warm_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_wet_quarter$X1)
  mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X1 %in% test,]
  colnames(mean_temp_wet_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_wet_quarter <- mean_temp_wet_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, min_temp_coldest_month$X1)
  min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X1 %in% test,]
  colnames(min_temp_coldest_month) <- c("orthogroup", "min_sdP_DS")
  min_temp_coldest_month <- min_temp_coldest_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, prec_clim_change$X1)
  prec_clim_change <- prec_clim_change[prec_clim_change$X1 %in% test,]
  colnames(prec_clim_change) <- c("orthogroup", "min_sdP_DS")
  prec_clim_change <- prec_clim_change %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_cold_quarter$X1)
  precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X1 %in% test,]
  colnames(precip_cold_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_cold_quarter <- precip_cold_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_dry_month$X1)
  precip_dry_month <- precip_dry_month[precip_dry_month$X1 %in% test,]
  colnames(precip_dry_month) <- c("orthogroup", "min_sdP_DS")
  precip_dry_month <- precip_dry_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_dry_quarter$X1)
  precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X1 %in% test,]
  colnames(precip_dry_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_dry_quarter <- precip_dry_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_seasonality$X1)
  precip_seasonality <- precip_seasonality[precip_seasonality$X1 %in% test,]
  colnames(precip_seasonality) <- c("orthogroup", "min_sdP_DS")
  precip_seasonality <- precip_seasonality %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, precip_warm_quarter$X1)
  precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X1 %in% test,]
  colnames(precip_warm_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_warm_quarter <- precip_warm_quarter %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, precip_wet_month$X1)
  precip_wet_month <- precip_wet_month[precip_wet_month$X1 %in% test,]
  colnames(precip_wet_month) <- c("orthogroup", "min_sdP_DS")
  precip_wet_month <- precip_wet_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_wet_quarter$X1)
  precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X1 %in% test,]
  colnames(precip_wet_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_wet_quarter <- precip_wet_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, temp_range$X1)
  temp_range <- temp_range[temp_range$X1 %in% test,]
  colnames(temp_range) <- c("orthogroup", "min_sdP_DS")
  temp_range <- temp_range %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, temp_seasonality$X1)
  temp_seasonality <- temp_seasonality[temp_seasonality$X1 %in% test,]
  colnames(temp_seasonality) <- c("orthogroup", "min_sdP_DS")
  temp_seasonality <- temp_seasonality %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, tmax_clim_change$X1)
  tmax_clim_change <- tmax_clim_change[tmax_clim_change$X1 %in% test,]
  colnames(tmax_clim_change) <- c("orthogroup", "min_sdP_DS")
  tmax_clim_change <- tmax_clim_change %>%
  left_join(gab, by = c("orthogroup"))


my_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)

pdf("kubota_ahalleri_scatter.pdf", height = 14, width = 14)
par(mfrow = c(7, 3), mar=c(4,4,4,4), cex.axis=0.7, cex.lab=1, cex.main=1, cex.sub=1,mgp=c(2,1,0) )
lapply(my_list, function(x){
    return(plot(-log10(as.numeric(x$min_sdP_DS)),-log10(x$dunnsidak_orthopvalues), xlab="Jim_WZA_ortho_p", ylab="Gabriele_SweeD_ortho_p"))
})

dev.off()



##### Preparing Gabriele orthogroup p values
##### murray_esid
gab <- read.table("murray_esid_final_orthogroup.txt", header=TRUE)




##### Preparing Jim ortho p values
jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Eucalyptus sideroxylon",]

jim <-jim[order(jim$Orthogroup),]

jim_ranks <- split(jim$min_sdP_DS, jim$climate)

jim_ranks <- data.frame(jim_ranks)




 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))




  test <- intersect(gab$orthogroup, annual_precip$X1)
  annual_precip <- annual_precip[annual_precip$X1 %in% test,]
  colnames(annual_precip) <- c("orthogroup", "min_sdP_DS")
  annual_precip <- annual_precip %>%
  left_join(gab, by = c("orthogroup"))

  test <- intersect(gab$orthogroup, isothermality$X1)
  isothermality <- isothermality[isothermality$X1 %in% test,]
  colnames(isothermality) <- c("orthogroup", "min_sdP_DS")
  isothermality <- isothermality %>%
  left_join(gab, by = c("orthogroup"))

  test <- intersect(gab$orthogroup, max_temp_warmest_month$X1)
  max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X1 %in% test,]
  colnames(max_temp_warmest_month) <- c("orthogroup", "min_sdP_DS")
  max_temp_warmest_month <- max_temp_warmest_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_diurnal$X1)
  mean_diurnal <- mean_diurnal[mean_diurnal$X1 %in% test,]
  colnames(mean_diurnal) <- c("orthogroup", "min_sdP_DS")
  mean_diurnal <- mean_diurnal %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp$X1)
  mean_temp <- mean_temp[mean_temp$X1 %in% test,]
  colnames(mean_temp) <- c("orthogroup", "min_sdP_DS")
  mean_temp <- mean_temp %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_cold_quarter$X1)
  mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X1 %in% test,]
  colnames(mean_temp_cold_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_cold_quarter <- mean_temp_cold_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_dry_quarter$X1)
  mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X1 %in% test,]
  colnames(mean_temp_dry_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_dry_quarter <- mean_temp_dry_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_warm_quarter$X1)
  mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X1 %in% test,]
  colnames(mean_temp_warm_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_warm_quarter <- mean_temp_warm_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_wet_quarter$X1)
  mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X1 %in% test,]
  colnames(mean_temp_wet_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_wet_quarter <- mean_temp_wet_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, min_temp_coldest_month$X1)
  min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X1 %in% test,]
  colnames(min_temp_coldest_month) <- c("orthogroup", "min_sdP_DS")
  min_temp_coldest_month <- min_temp_coldest_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, prec_clim_change$X1)
  prec_clim_change <- prec_clim_change[prec_clim_change$X1 %in% test,]
  colnames(prec_clim_change) <- c("orthogroup", "min_sdP_DS")
  prec_clim_change <- prec_clim_change %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_cold_quarter$X1)
  precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X1 %in% test,]
  colnames(precip_cold_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_cold_quarter <- precip_cold_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_dry_month$X1)
  precip_dry_month <- precip_dry_month[precip_dry_month$X1 %in% test,]
  colnames(precip_dry_month) <- c("orthogroup", "min_sdP_DS")
  precip_dry_month <- precip_dry_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_dry_quarter$X1)
  precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X1 %in% test,]
  colnames(precip_dry_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_dry_quarter <- precip_dry_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_seasonality$X1)
  precip_seasonality <- precip_seasonality[precip_seasonality$X1 %in% test,]
  colnames(precip_seasonality) <- c("orthogroup", "min_sdP_DS")
  precip_seasonality <- precip_seasonality %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, precip_warm_quarter$X1)
  precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X1 %in% test,]
  colnames(precip_warm_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_warm_quarter <- precip_warm_quarter %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, precip_wet_month$X1)
  precip_wet_month <- precip_wet_month[precip_wet_month$X1 %in% test,]
  colnames(precip_wet_month) <- c("orthogroup", "min_sdP_DS")
  precip_wet_month <- precip_wet_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_wet_quarter$X1)
  precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X1 %in% test,]
  colnames(precip_wet_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_wet_quarter <- precip_wet_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, temp_range$X1)
  temp_range <- temp_range[temp_range$X1 %in% test,]
  colnames(temp_range) <- c("orthogroup", "min_sdP_DS")
  temp_range <- temp_range %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, temp_seasonality$X1)
  temp_seasonality <- temp_seasonality[temp_seasonality$X1 %in% test,]
  colnames(temp_seasonality) <- c("orthogroup", "min_sdP_DS")
  temp_seasonality <- temp_seasonality %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, tmax_clim_change$X1)
  tmax_clim_change <- tmax_clim_change[tmax_clim_change$X1 %in% test,]
  colnames(tmax_clim_change) <- c("orthogroup", "min_sdP_DS")
  tmax_clim_change <- tmax_clim_change %>%
  left_join(gab, by = c("orthogroup"))


my_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)

pdf("murray_esid_scatter.pdf", height = 14, width = 14)
par(mfrow = c(7, 3), mar=c(4,4,4,4), cex.axis=0.7, cex.lab=1, cex.main=1, cex.sub=1,mgp=c(2,1,0) )
lapply(my_list, function(x){
    return(plot(-log10(as.numeric(x$min_sdP_DS)),-log10(x$dunnsidak_orthopvalues), xlab="Jim_WZA_ortho_p", ylab="Gabriele_SweeD_ortho_p"))
})

dev.off()



##### Preparing Gabriele orthogroup p values
##### savolainen_alyrata
gab <- read.table("savolainen_alyrata_final_orthogroup.txt", header=TRUE)




##### Preparing Jim ortho p values
jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Arabidopsis lyrata",]

jim <-jim[order(jim$Orthogroup),]

jim_ranks <- split(jim$min_sdP_DS, jim$climate)

jim_ranks <- data.frame(jim_ranks)




 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))




  test <- intersect(gab$orthogroup, annual_precip$X1)
  annual_precip <- annual_precip[annual_precip$X1 %in% test,]
  colnames(annual_precip) <- c("orthogroup", "min_sdP_DS")
  annual_precip <- annual_precip %>%
  left_join(gab, by = c("orthogroup"))

  test <- intersect(gab$orthogroup, isothermality$X1)
  isothermality <- isothermality[isothermality$X1 %in% test,]
  colnames(isothermality) <- c("orthogroup", "min_sdP_DS")
  isothermality <- isothermality %>%
  left_join(gab, by = c("orthogroup"))

  test <- intersect(gab$orthogroup, max_temp_warmest_month$X1)
  max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X1 %in% test,]
  colnames(max_temp_warmest_month) <- c("orthogroup", "min_sdP_DS")
  max_temp_warmest_month <- max_temp_warmest_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_diurnal$X1)
  mean_diurnal <- mean_diurnal[mean_diurnal$X1 %in% test,]
  colnames(mean_diurnal) <- c("orthogroup", "min_sdP_DS")
  mean_diurnal <- mean_diurnal %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp$X1)
  mean_temp <- mean_temp[mean_temp$X1 %in% test,]
  colnames(mean_temp) <- c("orthogroup", "min_sdP_DS")
  mean_temp <- mean_temp %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_cold_quarter$X1)
  mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X1 %in% test,]
  colnames(mean_temp_cold_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_cold_quarter <- mean_temp_cold_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_dry_quarter$X1)
  mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X1 %in% test,]
  colnames(mean_temp_dry_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_dry_quarter <- mean_temp_dry_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_warm_quarter$X1)
  mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X1 %in% test,]
  colnames(mean_temp_warm_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_warm_quarter <- mean_temp_warm_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_wet_quarter$X1)
  mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X1 %in% test,]
  colnames(mean_temp_wet_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_wet_quarter <- mean_temp_wet_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, min_temp_coldest_month$X1)
  min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X1 %in% test,]
  colnames(min_temp_coldest_month) <- c("orthogroup", "min_sdP_DS")
  min_temp_coldest_month <- min_temp_coldest_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, prec_clim_change$X1)
  prec_clim_change <- prec_clim_change[prec_clim_change$X1 %in% test,]
  colnames(prec_clim_change) <- c("orthogroup", "min_sdP_DS")
  prec_clim_change <- prec_clim_change %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_cold_quarter$X1)
  precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X1 %in% test,]
  colnames(precip_cold_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_cold_quarter <- precip_cold_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_dry_month$X1)
  precip_dry_month <- precip_dry_month[precip_dry_month$X1 %in% test,]
  colnames(precip_dry_month) <- c("orthogroup", "min_sdP_DS")
  precip_dry_month <- precip_dry_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_dry_quarter$X1)
  precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X1 %in% test,]
  colnames(precip_dry_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_dry_quarter <- precip_dry_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_seasonality$X1)
  precip_seasonality <- precip_seasonality[precip_seasonality$X1 %in% test,]
  colnames(precip_seasonality) <- c("orthogroup", "min_sdP_DS")
  precip_seasonality <- precip_seasonality %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, precip_warm_quarter$X1)
  precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X1 %in% test,]
  colnames(precip_warm_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_warm_quarter <- precip_warm_quarter %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, precip_wet_month$X1)
  precip_wet_month <- precip_wet_month[precip_wet_month$X1 %in% test,]
  colnames(precip_wet_month) <- c("orthogroup", "min_sdP_DS")
  precip_wet_month <- precip_wet_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_wet_quarter$X1)
  precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X1 %in% test,]
  colnames(precip_wet_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_wet_quarter <- precip_wet_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, temp_range$X1)
  temp_range <- temp_range[temp_range$X1 %in% test,]
  colnames(temp_range) <- c("orthogroup", "min_sdP_DS")
  temp_range <- temp_range %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, temp_seasonality$X1)
  temp_seasonality <- temp_seasonality[temp_seasonality$X1 %in% test,]
  colnames(temp_seasonality) <- c("orthogroup", "min_sdP_DS")
  temp_seasonality <- temp_seasonality %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, tmax_clim_change$X1)
  tmax_clim_change <- tmax_clim_change[tmax_clim_change$X1 %in% test,]
  colnames(tmax_clim_change) <- c("orthogroup", "min_sdP_DS")
  tmax_clim_change <- tmax_clim_change %>%
  left_join(gab, by = c("orthogroup"))


my_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)

pdf("savolainen_alyrata_scatter.pdf", height = 14, width = 14)
par(mfrow = c(7, 3), mar=c(4,4,4,4), cex.axis=0.7, cex.lab=1, cex.main=1, cex.sub=1,mgp=c(2,1,0) )
lapply(my_list, function(x){
    return(plot(-log10(as.numeric(x$min_sdP_DS)),-log10(x$dunnsidak_orthopvalues), xlab="Jim_WZA_ortho_p", ylab="Gabriele_SweeD_ortho_p"))
})

dev.off()



##### Preparing Gabriele orthogroup p values
##### mitchell_bstricta
gab <- read.table("mitchell_bstricta_final_orthogroup.txt", header=TRUE)




##### Preparing Jim ortho p values
jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Boechera stricta",]

jim <-jim[order(jim$Orthogroup),]

jim_ranks <- split(jim$min_sdP_DS, jim$climate)

jim_ranks <- data.frame(jim_ranks)




 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))




  test <- intersect(gab$orthogroup, annual_precip$X1)
  annual_precip <- annual_precip[annual_precip$X1 %in% test,]
  colnames(annual_precip) <- c("orthogroup", "min_sdP_DS")
  annual_precip <- annual_precip %>%
  left_join(gab, by = c("orthogroup"))

  test <- intersect(gab$orthogroup, isothermality$X1)
  isothermality <- isothermality[isothermality$X1 %in% test,]
  colnames(isothermality) <- c("orthogroup", "min_sdP_DS")
  isothermality <- isothermality %>%
  left_join(gab, by = c("orthogroup"))

  test <- intersect(gab$orthogroup, max_temp_warmest_month$X1)
  max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X1 %in% test,]
  colnames(max_temp_warmest_month) <- c("orthogroup", "min_sdP_DS")
  max_temp_warmest_month <- max_temp_warmest_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_diurnal$X1)
  mean_diurnal <- mean_diurnal[mean_diurnal$X1 %in% test,]
  colnames(mean_diurnal) <- c("orthogroup", "min_sdP_DS")
  mean_diurnal <- mean_diurnal %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp$X1)
  mean_temp <- mean_temp[mean_temp$X1 %in% test,]
  colnames(mean_temp) <- c("orthogroup", "min_sdP_DS")
  mean_temp <- mean_temp %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_cold_quarter$X1)
  mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X1 %in% test,]
  colnames(mean_temp_cold_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_cold_quarter <- mean_temp_cold_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_dry_quarter$X1)
  mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X1 %in% test,]
  colnames(mean_temp_dry_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_dry_quarter <- mean_temp_dry_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_warm_quarter$X1)
  mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X1 %in% test,]
  colnames(mean_temp_warm_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_warm_quarter <- mean_temp_warm_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_wet_quarter$X1)
  mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X1 %in% test,]
  colnames(mean_temp_wet_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_wet_quarter <- mean_temp_wet_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, min_temp_coldest_month$X1)
  min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X1 %in% test,]
  colnames(min_temp_coldest_month) <- c("orthogroup", "min_sdP_DS")
  min_temp_coldest_month <- min_temp_coldest_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, prec_clim_change$X1)
  prec_clim_change <- prec_clim_change[prec_clim_change$X1 %in% test,]
  colnames(prec_clim_change) <- c("orthogroup", "min_sdP_DS")
  prec_clim_change <- prec_clim_change %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_cold_quarter$X1)
  precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X1 %in% test,]
  colnames(precip_cold_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_cold_quarter <- precip_cold_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_dry_month$X1)
  precip_dry_month <- precip_dry_month[precip_dry_month$X1 %in% test,]
  colnames(precip_dry_month) <- c("orthogroup", "min_sdP_DS")
  precip_dry_month <- precip_dry_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_dry_quarter$X1)
  precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X1 %in% test,]
  colnames(precip_dry_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_dry_quarter <- precip_dry_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_seasonality$X1)
  precip_seasonality <- precip_seasonality[precip_seasonality$X1 %in% test,]
  colnames(precip_seasonality) <- c("orthogroup", "min_sdP_DS")
  precip_seasonality <- precip_seasonality %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, precip_warm_quarter$X1)
  precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X1 %in% test,]
  colnames(precip_warm_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_warm_quarter <- precip_warm_quarter %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, precip_wet_month$X1)
  precip_wet_month <- precip_wet_month[precip_wet_month$X1 %in% test,]
  colnames(precip_wet_month) <- c("orthogroup", "min_sdP_DS")
  precip_wet_month <- precip_wet_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_wet_quarter$X1)
  precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X1 %in% test,]
  colnames(precip_wet_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_wet_quarter <- precip_wet_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, temp_range$X1)
  temp_range <- temp_range[temp_range$X1 %in% test,]
  colnames(temp_range) <- c("orthogroup", "min_sdP_DS")
  temp_range <- temp_range %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, temp_seasonality$X1)
  temp_seasonality <- temp_seasonality[temp_seasonality$X1 %in% test,]
  colnames(temp_seasonality) <- c("orthogroup", "min_sdP_DS")
  temp_seasonality <- temp_seasonality %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, tmax_clim_change$X1)
  tmax_clim_change <- tmax_clim_change[tmax_clim_change$X1 %in% test,]
  colnames(tmax_clim_change) <- c("orthogroup", "min_sdP_DS")
  tmax_clim_change <- tmax_clim_change %>%
  left_join(gab, by = c("orthogroup"))


my_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)

pdf("mitchell_bstricta_scatter.pdf", height = 14, width = 14)
par(mfrow = c(7, 3), mar=c(4,4,4,4), cex.axis=0.7, cex.lab=1, cex.main=1, cex.sub=1,mgp=c(2,1,0) )
lapply(my_list, function(x){
    return(plot(-log10(as.numeric(x$min_sdP_DS)),-log10(x$dunnsidak_orthopvalues), xlab="Jim_WZA_ortho_p", ylab="Gabriele_SweeD_ortho_p"))
})

dev.off()


##### Preparing Gabriele orthogroup p values
##### murray_ealb
gab <- read.table("murray_ealb_final_orthogroup.txt", header=TRUE)




##### Preparing Jim ortho p values
jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Eucalyptus albens",]

jim <-jim[order(jim$Orthogroup),]

jim_ranks <- split(jim$min_sdP_DS, jim$climate)

jim_ranks <- data.frame(jim_ranks)




 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))




  test <- intersect(gab$orthogroup, annual_precip$X1)
  annual_precip <- annual_precip[annual_precip$X1 %in% test,]
  colnames(annual_precip) <- c("orthogroup", "min_sdP_DS")
  annual_precip <- annual_precip %>%
  left_join(gab, by = c("orthogroup"))

  test <- intersect(gab$orthogroup, isothermality$X1)
  isothermality <- isothermality[isothermality$X1 %in% test,]
  colnames(isothermality) <- c("orthogroup", "min_sdP_DS")
  isothermality <- isothermality %>%
  left_join(gab, by = c("orthogroup"))

  test <- intersect(gab$orthogroup, max_temp_warmest_month$X1)
  max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X1 %in% test,]
  colnames(max_temp_warmest_month) <- c("orthogroup", "min_sdP_DS")
  max_temp_warmest_month <- max_temp_warmest_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_diurnal$X1)
  mean_diurnal <- mean_diurnal[mean_diurnal$X1 %in% test,]
  colnames(mean_diurnal) <- c("orthogroup", "min_sdP_DS")
  mean_diurnal <- mean_diurnal %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp$X1)
  mean_temp <- mean_temp[mean_temp$X1 %in% test,]
  colnames(mean_temp) <- c("orthogroup", "min_sdP_DS")
  mean_temp <- mean_temp %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_cold_quarter$X1)
  mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X1 %in% test,]
  colnames(mean_temp_cold_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_cold_quarter <- mean_temp_cold_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_dry_quarter$X1)
  mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X1 %in% test,]
  colnames(mean_temp_dry_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_dry_quarter <- mean_temp_dry_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_warm_quarter$X1)
  mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X1 %in% test,]
  colnames(mean_temp_warm_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_warm_quarter <- mean_temp_warm_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_wet_quarter$X1)
  mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X1 %in% test,]
  colnames(mean_temp_wet_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_wet_quarter <- mean_temp_wet_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, min_temp_coldest_month$X1)
  min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X1 %in% test,]
  colnames(min_temp_coldest_month) <- c("orthogroup", "min_sdP_DS")
  min_temp_coldest_month <- min_temp_coldest_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, prec_clim_change$X1)
  prec_clim_change <- prec_clim_change[prec_clim_change$X1 %in% test,]
  colnames(prec_clim_change) <- c("orthogroup", "min_sdP_DS")
  prec_clim_change <- prec_clim_change %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_cold_quarter$X1)
  precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X1 %in% test,]
  colnames(precip_cold_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_cold_quarter <- precip_cold_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_dry_month$X1)
  precip_dry_month <- precip_dry_month[precip_dry_month$X1 %in% test,]
  colnames(precip_dry_month) <- c("orthogroup", "min_sdP_DS")
  precip_dry_month <- precip_dry_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_dry_quarter$X1)
  precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X1 %in% test,]
  colnames(precip_dry_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_dry_quarter <- precip_dry_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_seasonality$X1)
  precip_seasonality <- precip_seasonality[precip_seasonality$X1 %in% test,]
  colnames(precip_seasonality) <- c("orthogroup", "min_sdP_DS")
  precip_seasonality <- precip_seasonality %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, precip_warm_quarter$X1)
  precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X1 %in% test,]
  colnames(precip_warm_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_warm_quarter <- precip_warm_quarter %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, precip_wet_month$X1)
  precip_wet_month <- precip_wet_month[precip_wet_month$X1 %in% test,]
  colnames(precip_wet_month) <- c("orthogroup", "min_sdP_DS")
  precip_wet_month <- precip_wet_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_wet_quarter$X1)
  precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X1 %in% test,]
  colnames(precip_wet_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_wet_quarter <- precip_wet_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, temp_range$X1)
  temp_range <- temp_range[temp_range$X1 %in% test,]
  colnames(temp_range) <- c("orthogroup", "min_sdP_DS")
  temp_range <- temp_range %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, temp_seasonality$X1)
  temp_seasonality <- temp_seasonality[temp_seasonality$X1 %in% test,]
  colnames(temp_seasonality) <- c("orthogroup", "min_sdP_DS")
  temp_seasonality <- temp_seasonality %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, tmax_clim_change$X1)
  tmax_clim_change <- tmax_clim_change[tmax_clim_change$X1 %in% test,]
  colnames(tmax_clim_change) <- c("orthogroup", "min_sdP_DS")
  tmax_clim_change <- tmax_clim_change %>%
  left_join(gab, by = c("orthogroup"))


my_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)

pdf("murray_ealb_scatter.pdf", height = 14, width = 14)
par(mfrow = c(7, 3), mar=c(4,4,4,4), cex.axis=0.7, cex.lab=1, cex.main=1, cex.sub=1,mgp=c(2,1,0) )
lapply(my_list, function(x){
    return(plot(-log10(as.numeric(x$min_sdP_DS)),-log10(x$dunnsidak_orthopvalues), xlab="Jim_WZA_ortho_p", ylab="Gabriele_SweeD_ortho_p"))
})

dev.off()


##### Preparing Gabriele orthogroup p values
##### tiffin_truncatula
gab <- read.table("tiffin_truncatula_final_orthogroup.txt", header=TRUE)




##### Preparing Jim ortho p values
jim <- readRDS("run220927_COMBINED25species_mirrored_221026_WZA_OG_pvals.rds")
jim <- jim[jim$species == "Medicago truncatula",]

jim <-jim[order(jim$Orthogroup),]

jim_ranks <- split(jim$min_sdP_DS, jim$climate)

jim_ranks <- data.frame(jim_ranks)




 annual_precip <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$annual_precip))
 isothermality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$isothermality))
 max_temp_warmest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$max_temp_warmest_month))
 mean_diurnal <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_diurnal))
 mean_temp <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp))
 mean_temp_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_cold_quarter))
 mean_temp_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_dry_quarter))
 mean_temp_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_warm_quarter))
 mean_temp_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$mean_temp_wet_quarter))
 min_temp_coldest_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$min_temp_coldest_month))
 prec_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$prec_clim_change))
 precip_cold_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_cold_quarter))
 precip_dry_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_month))
 precip_dry_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_dry_quarter))
 precip_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_seasonality))
 precip_warm_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_warm_quarter))
 precip_wet_month <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_month))
 precip_wet_quarter <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$precip_wet_quarter))
 temp_range <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_range))
 temp_seasonality <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$temp_seasonality))
 tmax_clim_change <- data.frame(cbind(unique(jim$Orthogroup),jim_ranks$tmax_clim_change))




  test <- intersect(gab$orthogroup, annual_precip$X1)
  annual_precip <- annual_precip[annual_precip$X1 %in% test,]
  colnames(annual_precip) <- c("orthogroup", "min_sdP_DS")
  annual_precip <- annual_precip %>%
  left_join(gab, by = c("orthogroup"))

  test <- intersect(gab$orthogroup, isothermality$X1)
  isothermality <- isothermality[isothermality$X1 %in% test,]
  colnames(isothermality) <- c("orthogroup", "min_sdP_DS")
  isothermality <- isothermality %>%
  left_join(gab, by = c("orthogroup"))

  test <- intersect(gab$orthogroup, max_temp_warmest_month$X1)
  max_temp_warmest_month <- max_temp_warmest_month[max_temp_warmest_month$X1 %in% test,]
  colnames(max_temp_warmest_month) <- c("orthogroup", "min_sdP_DS")
  max_temp_warmest_month <- max_temp_warmest_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_diurnal$X1)
  mean_diurnal <- mean_diurnal[mean_diurnal$X1 %in% test,]
  colnames(mean_diurnal) <- c("orthogroup", "min_sdP_DS")
  mean_diurnal <- mean_diurnal %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp$X1)
  mean_temp <- mean_temp[mean_temp$X1 %in% test,]
  colnames(mean_temp) <- c("orthogroup", "min_sdP_DS")
  mean_temp <- mean_temp %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_cold_quarter$X1)
  mean_temp_cold_quarter <- mean_temp_cold_quarter[mean_temp_cold_quarter$X1 %in% test,]
  colnames(mean_temp_cold_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_cold_quarter <- mean_temp_cold_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_dry_quarter$X1)
  mean_temp_dry_quarter <- mean_temp_dry_quarter[mean_temp_dry_quarter$X1 %in% test,]
  colnames(mean_temp_dry_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_dry_quarter <- mean_temp_dry_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_warm_quarter$X1)
  mean_temp_warm_quarter <- mean_temp_warm_quarter[mean_temp_warm_quarter$X1 %in% test,]
  colnames(mean_temp_warm_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_warm_quarter <- mean_temp_warm_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, mean_temp_wet_quarter$X1)
  mean_temp_wet_quarter <- mean_temp_wet_quarter[mean_temp_wet_quarter$X1 %in% test,]
  colnames(mean_temp_wet_quarter) <- c("orthogroup", "min_sdP_DS")
  mean_temp_wet_quarter <- mean_temp_wet_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, min_temp_coldest_month$X1)
  min_temp_coldest_month <- min_temp_coldest_month[min_temp_coldest_month$X1 %in% test,]
  colnames(min_temp_coldest_month) <- c("orthogroup", "min_sdP_DS")
  min_temp_coldest_month <- min_temp_coldest_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, prec_clim_change$X1)
  prec_clim_change <- prec_clim_change[prec_clim_change$X1 %in% test,]
  colnames(prec_clim_change) <- c("orthogroup", "min_sdP_DS")
  prec_clim_change <- prec_clim_change %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_cold_quarter$X1)
  precip_cold_quarter <- precip_cold_quarter[precip_cold_quarter$X1 %in% test,]
  colnames(precip_cold_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_cold_quarter <- precip_cold_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_dry_month$X1)
  precip_dry_month <- precip_dry_month[precip_dry_month$X1 %in% test,]
  colnames(precip_dry_month) <- c("orthogroup", "min_sdP_DS")
  precip_dry_month <- precip_dry_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_dry_quarter$X1)
  precip_dry_quarter <- precip_dry_quarter[precip_dry_quarter$X1 %in% test,]
  colnames(precip_dry_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_dry_quarter <- precip_dry_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_seasonality$X1)
  precip_seasonality <- precip_seasonality[precip_seasonality$X1 %in% test,]
  colnames(precip_seasonality) <- c("orthogroup", "min_sdP_DS")
  precip_seasonality <- precip_seasonality %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, precip_warm_quarter$X1)
  precip_warm_quarter <- precip_warm_quarter[precip_warm_quarter$X1 %in% test,]
  colnames(precip_warm_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_warm_quarter <- precip_warm_quarter %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, precip_wet_month$X1)
  precip_wet_month <- precip_wet_month[precip_wet_month$X1 %in% test,]
  colnames(precip_wet_month) <- c("orthogroup", "min_sdP_DS")
  precip_wet_month <- precip_wet_month %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, precip_wet_quarter$X1)
  precip_wet_quarter <- precip_wet_quarter[precip_wet_quarter$X1 %in% test,]
  colnames(precip_wet_quarter) <- c("orthogroup", "min_sdP_DS")
  precip_wet_quarter <- precip_wet_quarter %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, temp_range$X1)
  temp_range <- temp_range[temp_range$X1 %in% test,]
  colnames(temp_range) <- c("orthogroup", "min_sdP_DS")
  temp_range <- temp_range %>%
  left_join(gab, by = c("orthogroup"))



  test <- intersect(gab$orthogroup, temp_seasonality$X1)
  temp_seasonality <- temp_seasonality[temp_seasonality$X1 %in% test,]
  colnames(temp_seasonality) <- c("orthogroup", "min_sdP_DS")
  temp_seasonality <- temp_seasonality %>%
  left_join(gab, by = c("orthogroup"))


  test <- intersect(gab$orthogroup, tmax_clim_change$X1)
  tmax_clim_change <- tmax_clim_change[tmax_clim_change$X1 %in% test,]
  colnames(tmax_clim_change) <- c("orthogroup", "min_sdP_DS")
  tmax_clim_change <- tmax_clim_change %>%
  left_join(gab, by = c("orthogroup"))


my_list <- list(annual_precip,isothermality, max_temp_warmest_month, mean_diurnal, mean_temp, mean_temp_cold_quarter, mean_temp_dry_quarter, mean_temp_warm_quarter, mean_temp_wet_quarter, min_temp_coldest_month, prec_clim_change,precip_cold_quarter,precip_dry_month,precip_dry_quarter,precip_seasonality,precip_warm_quarter,precip_wet_month,precip_wet_quarter,temp_range,temp_seasonality,tmax_clim_change)

pdf("tiffin_truncatula_scatter.pdf", height = 14, width = 14)
par(mfrow = c(7, 3), mar=c(4,4,4,4), cex.axis=0.7, cex.lab=1, cex.main=1, cex.sub=1,mgp=c(2,1,0) )
lapply(my_list, function(x){
    return(plot(-log10(as.numeric(x$min_sdP_DS)),-log10(x$dunnsidak_orthopvalues), xlab="Jim_WZA_ortho_p", ylab="Gabriele_SweeD_ortho_p"))
})

dev.off()
