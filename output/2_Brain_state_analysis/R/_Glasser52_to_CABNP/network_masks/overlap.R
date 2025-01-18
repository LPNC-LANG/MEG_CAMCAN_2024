################################################################################
# Written by Clément Guichet, PhD Student
# LPNC - CNRS UMR 5105
# 2024

################################################################################

library(tidyverse)
library(data.table)
library(plyr)
library(janitor)
library(jsonlite)
library(stringr)

rm(list = ls())
setwd("E:/Research_Projects/MEG_CamCAN/network_masks/glasser52/")

# With the full 360-region HCP atlas ----
region_list <- rio::import("region_list52.txt")
region_list_clean <- rbind(colnames(region_list), region_list)
colnames(region_list_clean) <- "region"

library(readr)
tmp <- read_table("../output/MD_overlap.tsv", col_names = FALSE) %>% dplyr::select(X9)
MD_overlap <- tmp %>%
  mutate(MD_overlap = X9 / lag(X9) * 100) %>% 
  .[-seq(1, nrow(.), 2), 2] %>% 
  as.data.frame() %>%
  mutate(parcel_number = region_list_clean$region) %>% 
  dplyr::mutate(across("parcel_number", substr, 8, str_length(parcel_number))) %>% 
  arrange(MD_overlap, decreasing = T)

tmp <- read_table("../output/SCN_overlap.tsv", col_names = FALSE) %>% dplyr::select(X9)
SCN_overlap <- tmp %>%
  mutate(SCN_overlap = X9 / lag(X9) * 100) %>% 
  .[-seq(1, nrow(.), 2), 2] %>% 
  as.data.frame() %>%
  mutate(parcel_number = region_list_clean$region) %>% 
  dplyr::mutate(across("parcel_number", substr, 8, str_length(parcel_number))) %>% 
  arrange(SCN_overlap, decreasing = T)

tmp <- read_table("../output/DMN_overlap.tsv", col_names = FALSE) %>% dplyr::select(X9)
DMN_overlap <- tmp %>%
  mutate(DMN_overlap = X9 / lag(X9) * 100) %>% 
  .[-seq(1, nrow(.), 2), 2] %>% 
  as.data.frame() %>%
  mutate(parcel_number = region_list_clean$region) %>% 
  dplyr::mutate(across("parcel_number", substr, 8, str_length(parcel_number))) %>% 
  arrange(DMN_overlap, decreasing = T)

tmp <- read_table("../output/SMN_overlap.tsv", col_names = FALSE) %>% dplyr::select(X9)
SMN_overlap <- tmp %>%
  mutate(SMN_overlap = X9 / lag(X9) * 100) %>% 
  .[-seq(1, nrow(.), 2), 2] %>% 
  as.data.frame() %>%
  mutate(parcel_number = region_list_clean$region) %>% 
  dplyr::mutate(across("parcel_number", substr, 8, str_length(parcel_number))) %>% 
  arrange(SMN_overlap, decreasing = T)

tmp <- read_table("../output/LANG_overlap.tsv", col_names = FALSE) %>% dplyr::select(X9)
LANG_overlap <- tmp %>%
  mutate(LANG_overlap = X9 / lag(X9) * 100) %>% 
  .[-seq(1, nrow(.), 2), 2] %>% 
  as.data.frame() %>%
  mutate(parcel_number = region_list_clean$region) %>% 
  dplyr::mutate(across("parcel_number", substr, 8, str_length(parcel_number))) %>% 
  arrange(LANG_overlap, decreasing = T)

tmp <- read_table("../output/CON_overlap.tsv", col_names = FALSE) %>% dplyr::select(X9)
CON_overlap <- tmp %>%
  mutate(CON_overlap = X9 / lag(X9) * 100) %>% 
  .[-seq(1, nrow(.), 2), 2] %>% 
  as.data.frame() %>%
  mutate(parcel_number = region_list_clean$region) %>% 
  dplyr::mutate(across("parcel_number", substr, 8, str_length(parcel_number))) %>% 
  arrange(CON_overlap, decreasing = T)

tmp <- read_table("../output/FPN_overlap.tsv", col_names = FALSE) %>% dplyr::select(X9)
FPN_overlap <- tmp %>%
  mutate(FPN_overlap = X9 / lag(X9) * 100) %>% 
  .[-seq(1, nrow(.), 2), 2] %>% 
  as.data.frame() %>%
  mutate(parcel_number = region_list_clean$region) %>% 
  dplyr::mutate(across("parcel_number", substr, 8, str_length(parcel_number))) %>% 
  arrange(FPN_overlap, decreasing = T)

HCP_regridded <- rio::import("D:/Analyses/Atlases/HCP/LMN_HCP/1_Isolate_HCP52_ROIs/parcel_names_and_mni_coordinates.xlsx")
RSN_overlap <- merge(HCP_regridded, MD_overlap, by = "parcel_number") %>% 
  merge(., SCN_overlap, by = "parcel_number") %>% 
  merge(., DMN_overlap, by = "parcel_number") %>% 
  merge(., SMN_overlap, by = "parcel_number") %>% 
  merge(., LANG_overlap, by = "parcel_number") %>% 
  merge(., CON_overlap, by = "parcel_number") %>% 
  merge(., FPN_overlap, by = "parcel_number")


# Overlap with power maps -------------------------------------------------
library(RcppCNPy)
group_level_psd <- RcppCNPy::npyLoad("E:/Research_Projects/MEG_CamCAN/TDE_HMM/results/08_states_networks/grp_psd_centered.npy") %>% 
  t() %>% as.data.frame()
rownames(group_level_psd) <- HCP_regridded$Parcel

nparcels = 20
keep_top <- function(column) {
  # Get the indices
  top_indices <- order(column, decreasing = TRUE)[1:nparcels]
  result <- numeric(length(column))
  result[top_indices] <- column[top_indices]
  return(result)
}

group_level_psd_top <- as.data.frame(apply(group_level_psd, 2, keep_top))


RSN_by_parcel <- t(as.matrix(RSN_overlap[,6:12]))
parcel_by_state <- as.matrix(group_level_psd_top)

RSN_maps <- RSN_by_parcel %*% parcel_by_state %>% t() %>% scale(center=F) %>% t() %>% as.data.frame() 
# 
# plot_radar_RSN <- function(i){
#   library(tidyverse)
#   library(stringr)
#   library(scales)
#   library(prismatic)
#   
#   toplot <- RSN_maps[,i] %>% as.data.frame() %>% mutate(RSN = colnames(RSN_overlap)[6:17])
#   colnames(toplot) <- c("V1", "RSN")
#   
#   plot_alpha <- ggplot(toplot) +
#     
#     #make custom panel grid
#     geom_hline(yintercept = min(toplot$V1), color = "lightgrey") +
#     geom_hline(yintercept = quantile(toplot$V1)[2], color = "lightgrey") +
#     geom_hline(yintercept = 0, color = "gray12", size = 2, linetype = 4) +
#     geom_hline(yintercept = quantile(toplot$V1)[3], color = "lightgrey") +
#     geom_hline(yintercept = max(toplot$V1), color = "lightgrey") +
#     
#     geom_col(aes(
#       x = reorder(RSN,V1), 
#       y = V1,#is numeric
#       fill = V1), 
#       position = "dodge2",
#       show.legend = TRUE,
#       alpha = .9) +
#     
#     if (max(toplot$V1) <= 0){ # If only negative
#       scale_fill_distiller(palette = "Blues")
#     } else if (min(toplot$V1) >= 0){ # If only positive
#       scale_fill_distiller(palette = "Reds", direction = 1)
#     } else if (min(toplot$V1) <= 0 && max(toplot$V1) >= 0){ # If both
#       scale_fill_gradientn("Power (a.u.)",
#                            colours = c("#0571B0FF", "#92C5DEFF", "#F4A582FF", "#CA0020FF"))
#     }
#   
#   
#   plot_alpha + geom_point(aes(x = reorder(RSN,V1),
#                               y = V1),
#                           size = 3,
#                           color = "gray12")+
#     
#     #lollipop shaft
#     geom_segment(aes(
#       x = reorder(RSN,V1),
#       y = min(V1),
#       xend = reorder(RSN,V1),
#       yend = max(V1)),
#       linetype = "dashed",
#       color = "gray12") +
#     
#     #transform to polar coordinate system
#     coord_polar() +
#     
#     #theming
#     theme(legend.position = "none",
#           axis.title = element_blank(),
#           axis.ticks = element_blank(),
#           axis.text.y = element_blank(),
#           axis.text.x = element_text(color = "gray12",
#                                      size = 20),
#           panel.background = element_rect(fill = "white",
#                                           color = "white"),
#           panel.grid = element_blank(),
#           panel.grid.major.x = element_blank(),
#           text = element_text(color = "gray12",
#                               family = "Arial"))
#   
#   ggsave(path = r"(../output)",
#          filename = sprintf("State_%d_RSN_composition.png",i),
#          width = 20,
#          height = 20,
#          units = "cm",
#          type = "cairo-png",
#          dpi = 400)
# }
# 
# lapply(seq(8), plot_radar_RSN)
