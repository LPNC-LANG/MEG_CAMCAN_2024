# Script for formatting output from network overlap calculator
# CG - November 2024

################################################################################
################################################################################

library(tidyverse)
library(jsonlite)
library(plyr)
library(data.table)
library(janitor)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#Read file from JSON
overlap <- as.data.frame(fromJSON('./3_Multiplied/Overlap_voxels.json')) #JSON for the liberal CAB-NP


Label <- read.delim("region_list52.txt", header = F) %>% 
  plyr::rename(c("V1" = "parcel_number"))

# Data manipulation------------------------------------------------------------
#Change accordingly
version <- overlap

# Rbind every two columns------------------------------------------------------
lst <- split.default(version, cumsum(rep(c(TRUE, FALSE), ncol(overlap)/2))) 
data <- rbindlist(setNames(lst, seq_along(lst)), idcol="idx", use.names = FALSE) %>% 
  mutate_at("idx", funs(as.numeric(as.character(.)))) %>% # This is the alphabetical order (1,10,11,...) used by MATLAB
  plyr::rename(c("V1" = "Network")) %>% 
  plyr::rename(c("V2" = "voxel_proportion"))
  
data_wider <- data %>% 
  #Create a unique identifier row
  group_by(idx, Network) %>%
  dplyr::mutate(row = row_number()) %>%
  pivot_wider(
    names_from = "Network",
    values_from = "voxel_proportion"
  ) %>% 
  subset(row == "1") %>% 
  dplyr::select(-c(row, `0`)) %>% 
  replace(is.na(.), 0) %>% 
  
  plyr::rename(c("1" = "Visual 1")) %>% 
  plyr::rename(c("2" = "Visual 2")) %>% 
  plyr::rename(c("3" = "SMN")) %>% 
  plyr::rename(c("4" = "CON")) %>% 
  plyr::rename(c("5" = "DAN")) %>% 
  plyr::rename(c("6" = "Language")) %>% 
  plyr::rename(c("7" = "FPN")) %>%
  plyr::rename(c("8" = "Auditory")) %>% 
  plyr::rename(c("9" = "DMN")) %>% 
  plyr::rename(c("10" = "PMM")) %>% 
  plyr::rename(c("11" = "VMM")) %>% 
  plyr::rename(c("12" = "OA"))

###################################################################################
setwd("E:/Research_Projects/MEG_CamCAN/TDE_HMM/")

###################################################################################
HCP_regridded <- rio::import("./output/2_Brain_state_analysis/R/_Glasser52_to_CABNP/additional_info/parcel_names_and_mni_coordinates.xlsx")

RSN_overlap <- data_wider %>% cbind(Label, .) %>% .[,-2] %>% 
  dplyr::mutate(across("parcel_number", substr, 8, str_length(parcel_number))) %>% 
  mutate_at(vars(parcel_number), funs(as.numeric(.))) %>% 
  merge(HCP_regridded, ., by = "parcel_number")

library(writexl)
write_xlsx(RSN_overlap, "./output/2_Brain_state_analysis/R/_Glasser52_to_CABNP/output/Glasser52_RSN_composition.xlsx")

# PSD-weighted RSN composition of each brain state
library(RcppCNPy)
group_level_psd <- RcppCNPy::npyLoad("./output/2_Brain_state_analysis/R/gpow_dyn.npy") %>% 
  t() %>% as.data.frame()
rownames(group_level_psd) <- HCP_regridded$Parcel
write.csv(group_level_psd*1000, "./output/2_Brain_state_analysis/R/_Glasser52_to_CABNP/output/region_PSD_maps.csv")

nparcels = 52
keep_top <- function(column) {
  # Get the indices
  top_indices <- order(column, decreasing = TRUE)[1:nparcels]
  
  # Assign the top 10 values to their respective positions
  result <- numeric(length(column))
  result[top_indices] <- column[top_indices]
  return(result)
}
group_level_psd_top <- as.data.frame(apply(group_level_psd, 2, keep_top))


RSN_by_parcel <- t(as.matrix(RSN_overlap[,6:17]))
parcel_by_state <- as.matrix(group_level_psd_top)

RSN_maps <- RSN_by_parcel %*% parcel_by_state %>% as.data.frame() 
write.csv(RSN_maps, "./output/2_Brain_state_analysis/R/_Glasser52_to_CABNP/output/network_PSD_maps.csv")


plot_radar_RSN <- function(i){
  library(tidyverse)
  library(stringr)
  library(scales)
  library(prismatic)
  
  toplot <- RSN_maps[,i] %>% as.data.frame() %>% mutate(RSN = colnames(RSN_overlap)[6:17])
  colnames(toplot) <- c("V1", "RSN")
  
  plot_alpha <- ggplot(toplot) +
    
    #make custom panel grid
    geom_hline(yintercept = min(toplot$V1), color = "lightgrey") +
    geom_hline(yintercept = quantile(toplot$V1)[2], color = "lightgrey") +
    geom_hline(yintercept = 0, color = "gray12", size = 2, linetype = 4) +
    geom_hline(yintercept = quantile(toplot$V1)[3], color = "lightgrey") +
    geom_hline(yintercept = max(toplot$V1), color = "lightgrey") +
    
    geom_col(aes(
      x = reorder(RSN,V1), 
      y = V1,#is numeric
      fill = V1), 
      position = "dodge2",
      show.legend = TRUE,
      alpha = .9) +
  
    if (max(toplot$V1) <= 0){ # If only negative
      scale_fill_distiller(palette = "Blues")
    } else if (min(toplot$V1) >= 0){ # If only positive
        scale_fill_distiller(palette = "Reds", direction = 1)
    } else if (min(toplot$V1) <= 0 && max(toplot$V1) >= 0){ # If both
        scale_fill_gradientn("Power (a.u.)",
                           colours = c("#0571B0FF", "#92C5DEFF", "#F4A582FF", "#CA0020FF"))
    }
    
    
    plot_alpha + geom_point(aes(x = reorder(RSN,V1),
                   y = V1),
               size = 3,
               color = "gray12")+
    
    #lollipop shaft
    geom_segment(aes(
      x = reorder(RSN,V1),
      y = min(V1),
      xend = reorder(RSN,V1),
      yend = max(V1)),
      linetype = "dashed",
      color = "gray12") +
    
    #transform to polar coordinate system
    coord_polar() +
    
    #theming
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(color = "gray12",
                                     size = 20),
          panel.background = element_rect(fill = "white",
                                          color = "white"),
          panel.grid = element_blank(),
          panel.grid.major.x = element_blank(),
          text = element_text(color = "gray12",
                              family = "Arial"))
    
  ggsave(path = r"(./output/2_Brain_state_analysis/R/_Glasser52_to_CABNP/output/)",
         filename = sprintf("State_%d_RSN_composition.png",i),
         width = 20,
         height = 20,
         units = "cm",
         type = "cairo-png",
         dpi = 400)
}

lapply(seq(8), plot_radar_RSN)

  