
### CONTINUOUS HEATMAP
# library(ggplot2)
# library(reshape2)
# library(dplyr)
#
# df <- readRDS("orthogroup_results.rds")
# df <- df[, -c(1, 2, 3, 4, 6, 9)]
#
#
#
# # Define the custom order of the species
# custom_order <- c("ptricho", "ptremula", "bpendula", "bplaty","ealb", "emag", "esid", "mtruncatula", "athaliana", "ahalleri","capsella","bstricta","hann","hargo","hpet","atuber","phalli")
# custom_order <- rev(custom_order)
# # Convert the species column to a factor with the custom order
# df$species <- factor(df$species, levels = custom_order)
#
#
#
#   myplot<- ggplot(df, aes(x =reorder(Orthogroup, -table(Orthogroup)[Orthogroup]), y = species, fill = ortho_DS)) +
#   geom_tile() +
#   scale_fill_viridis_c(direction = -1, option  = "magma", breaks = c(0, 0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
#   theme_minimal() +
#   theme(axis.text.x = element_blank(),
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.direction = "vertical",
#         legend.position = "right",
#         legend.title = element_blank(),
#         axis.title.x = element_blank(),
#         axis.line.x = element_blank(),
#         axis.ticks.x = element_blank())
#
#
#
#
#
# ### DISCRETE HEATMAP
# library(ggplot2)
# library(reshape2)
# library(dplyr)
#
# df <- readRDS("orthogroup_results.rds")
# df <- df[, -c(1, 2, 3, 4, 6, 9)]
#
# df$ortho_DS <- ifelse(df$ortho_DS < 0.1, "sig", "not_sig")
# count_red <- df %>%   group_by(Orthogroup) %>%   summarise(count_red = sum(ortho_DS == "sig"))
# count_red <- as.data.frame(count_red)
# df <- df %>% left_join(count_red, by = c("Orthogroup"))
# df <- df[order(-df$count_red), ]
# df$Orthogroup <- factor(df$Orthogroup, levels = unique(df$Orthogroup))
# # Just plot OG with at least 1 ortho_DS < 0.1  -- "sig"
# df <- df[df$count_red > 0,]
#
#
# # Define the custom order of the species
# custom_order <- c("ptricho", "ptremula", "bpendula", "bplaty","ealb", "emag", "esid", "mtruncatula", "athaliana", "ahalleri","capsella","bstricta","hann","hargo","hpet","atuber","phalli")
# custom_order <- rev(custom_order)
# # Convert the species column to a factor with the custom order
# df$species <- factor(df$species, levels = custom_order)
#
#
# myplot2<- ggplot(df, aes(x =Orthogroup, y = species, fill = ortho_DS)) +
#   geom_tile() +
#   scale_fill_manual(values = c("navy", "gold")) +
#   theme_minimal() +
#   theme(axis.text.x = element_blank(),
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
# #       legend.direction = "vertical",
# #       legend.position = "right",
# #       legend.title = element_blank(),
#         legend.position = "none",
#         legend.title = element_blank(),
#         legend.key = element_blank(),
#         axis.title.x = element_blank(),
#         axis.line.x = element_blank(),
#         axis.ticks.x = element_blank())
#
#
#
#
# ### ONLY FDR < 0.4
# library(ggplot2)
# library(reshape2)
# library(dplyr)
#
# z <- readRDS("gab_picmin_results.rds")
# z <- z$picmin_res
# z <- z[z$picmin_fdr<0.4,]
# z <- as.data.frame(cbind(z$Orthogroup, z$climate_var))
# colnames(z) <- c("Orthogroup", "removeme")
#
#
# df <- readRDS("orthogroup_results.rds")
# df <- df %>%   left_join(z, by = c("Orthogroup"))
# df <- df[complete.cases(df),]
# df <- df[, -c(1, 2, 3, 4, 6, 9,10)]
#
# df$ortho_DS <- ifelse(df$ortho_DS < 0.1, "sig", "not_sig")
# count_red <- df %>%   group_by(Orthogroup) %>%   summarise(count_red = sum(ortho_DS == "sig"))
# count_red <- as.data.frame(count_red)
# df <- df %>% left_join(count_red, by = c("Orthogroup"))
# df <- df[order(-df$count_red), ]
# df$Orthogroup <- factor(df$Orthogroup, levels = unique(df$Orthogroup))
#
#
# # Define the custom order of the species
# custom_order <- c("ptricho", "ptremula", "bpendula", "bplaty","ealb", "emag", "esid", "mtruncatula", "athaliana", "ahalleri","capsella","bstricta","hann","hargo","hpet","atuber","phalli")
# custom_order <- rev(custom_order)
# # Convert the species column to a factor with the custom order
# df$species <- factor(df$species, levels = custom_order)
#
#
# myplot3<- ggplot(df, aes(x =Orthogroup, y = species, fill = ortho_DS)) +
#   geom_tile() +
#   scale_fill_manual(values = c("navy", "gold")) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
#   axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
# #       legend.direction = "vertical",
# #       legend.position = "right",
# #       legend.title = element_blank(),
#         legend.position = "none",
#         legend.title = element_blank(),
#         legend.key = element_blank(),
#         axis.title.x = element_blank(),
#         axis.line.x = element_blank(),
#         axis.ticks.x = element_blank())
#
#
#
# # tree
# library(cowplot)
# library(ggtree)
#
#
# tree <- read.tree("SpeciesTree_rooted_node_labels.txt")
# mytree_plot <- ggtree(tree)
#
# plot_grid(plot_grid(mytree_plot,ncol = 1), myplot, rel_widths = c(0.5, 1), nrow = 1, align = "h", axis = "tb")
# plot_grid(plot_grid(mytree_plot,ncol = 1), myplot2, rel_widths = c(0.5, 1), nrow = 1, align = "h", axis = "tb")
# plot_grid(plot_grid(mytree_plot,ncol = 1), myplot3, rel_widths = c(0.5, 1), nrow = 1, align = "h", axis = "tb")
#
#



### ONLY FDR < 0.4 using Heatmap
library(ggplot2)
library(reshape2)
library(dplyr)

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

heatmap.2(matrix_df, Rowv = NA, trace = "none", density.info = "none", col = Colors, breaks = Breaks, margins = c(10,10),keysize=1)

Colors=c("#9E0142","grey","grey","grey","grey","grey","grey","grey","grey","grey")
heatmap.2(matrix_df, Rowv = NA, trace = "none", density.info = "none", col = Colors, breaks = Breaks, margins = c(10,10), keysize=1)
