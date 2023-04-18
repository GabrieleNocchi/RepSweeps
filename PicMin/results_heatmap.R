library(ggplot2)
library(reshape2)
library(dplyr)
library(circlize)

z <- readRDS("gab_picmin_results.rds")
z <- z$picmin_res
z <- z[z$picmin_fdr<0.5,]
z <- as.data.frame(cbind(z$Orthogroup, z$climate_var))
colnames(z) <- c("Orthogroup", "removeme")


df <- readRDS("orthogroup_results.rds")
df <- df %>%   left_join(z, by = c("Orthogroup"))
df <- df[complete.cases(df),]
df <- df[, -c(2, 3, 4,5, 6, 9,10)]

df$species <- sub("ptricho", "P. trichocarpa", df$species)
df$species <- sub("emag", "E. magnificata", df$species)
df$species <- sub("ealb", "E. albens", df$species)
df$species <- sub("esid", "E. sideroxylon", df$species)
df$species <- sub("bplaty", "B. platyphylla", df$species)
df$species <- sub("bpendula", "B. pendula", df$species)
df$species <- sub("hann", "H. annus", df$species)
df$species <- sub("hargo", "H. argophyllus", df$species)
df$species <- sub("hpet", "H. petiolaris", df$species)
df$species <- sub("capsella", "C. rubella", df$species)
df$species <- sub("ahalleri", "A. halleri", df$species)
df$species <- sub("athaliana", "A. thaliana", df$species)
df$species <- sub("ptremula", "P. tremula", df$species)
df$species <- sub("bstricta", "B. stricta", df$species)
df$species <- sub("phalli", "P. hallii", df$species)
df$species <- sub("atuber", "A. tuberculatus", df$species)
df$species <- sub("mtruncatula", "M. truncatula", df$species)

# Define the custom order of the species
custom_order <- c("P. trichocarpa", "P. tremula", "B. pendula", "B. platyphylla","E. albens", "E. magnificata", "E. sideroxylon", "M. truncatula", "A. thaliana", "A. halleri","C. rubella","B. stricta","H. annus","H. argophyllus","H. petiolaris","A. tuberculatus","P. hallii")
custom_order <- rev(custom_order)
# Convert the species column to a factor with the custom order
df$species <- factor(df$species, levels = custom_order)
df$ortho_DS <- as.numeric(df$ortho_DS)

matrix_df <- dcast(df, species ~ Orthogroup, value.var = "ortho_DS")
rownames(matrix_df) <- matrix_df$species
matrix_df$species <- NULL
matrix_df <- as.matrix(matrix_df)



library(RColorBrewer)
library(gplots)
Colors=brewer.pal(11,"RdYlBu")
Colors=colorRampPalette(Colors)(10)
Breaks=seq(0,1,0.1)



library(ComplexHeatmap)
my_fdr <- as.data.frame(colnames(matrix_df))
colnames(my_fdr) <- "Orthogroup"
my_picmin <- readRDS("gab_picmin_results.rds")
my_picmin <- my_picmin$picmin_res
my_picmin <- my_picmin[, -c(2, 3, 5, 6)]
FDR <- my_fdr %>%   left_join(my_picmin, by = c("Orthogroup"))
column_ha = HeatmapAnnotation(`PicMin FDR` = anno_barplot(FDR$picmin_fdr, border = FALSE), annotation_name_side = "right", annotation_name_gp= gpar(fontsize = 9))
ht <- Heatmap(matrix_df, row_names_gp = gpar(fontsize = 8, fontface = "bold"), column_names_gp = gpar(fontsize = 10),col = Colors, na_col = "white", clustering_method_columns = "complete", cluster_rows = F, top_annotation = column_ha, show_column_dend = FALSE, heatmap_legend_param = list(title = "Empirical p"))
draw(ht, padding = unit(c(20, 20, 20, 20), "mm"))


#col_fun <- colorRamp2(c(0, 0.1, 1), c("#9E0142","grey","grey"))
ht <- Heatmap(matrix_df,row_names_gp = gpar(fontsize = 8, fontface = "bold"), column_names_gp = gpar(fontsize = 10), clustering_method_columns = "complete", col = c("#9E0142", "grey", "grey","grey","grey","grey" ,"grey", "grey","grey","grey"), na_col = "white", cluster_rows = F, top_annotation = column_ha, show_column_dend = FALSE, heatmap_legend_param = list(title = "Empirical p"))
draw(ht, padding = unit(c(20, 20, 20, 20), "mm" ))
