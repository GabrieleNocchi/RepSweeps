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



# Define the custom order of the species
custom_order <- c("ptricho", "ptremula", "bpendula", "bplaty","ealb", "emag", "esid", "mtruncatula", "athaliana", "ahalleri","capsella","bstricta","hann","hargo","hpet","atuber","phalli")
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
column_ha = HeatmapAnnotation(PicMin_FDR = anno_barplot(FDR$picmin_fdr, border = FALSE))
ht <- Heatmap(matrix_df, col = Colors, na_col = "white", cluster_rows = F, top_annotation = column_ha, show_column_dend = FALSE, heatmap_legend_param = list(title = "Emp-p"))
draw(ht, padding = unit(c(20, 20, 20, 20), "mm"))


col_fun <- colorRamp2(c(0, 0.1, 1), c("#9E0142", "grey", "grey"))
ht <- Heatmap(matrix_df, col = col_fun, na_col = "white", cluster_rows = F, top_annotation = column_ha, show_column_dend = FALSE, heatmap_legend_param = list(title = "Emp-p"))
draw(ht, padding = unit(c(20, 20, 20, 20), "mm"))
