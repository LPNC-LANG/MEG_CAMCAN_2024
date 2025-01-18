################################################################################
# Written by Clément Guichet, PhD Student
# LPNC - CNRS UMR 5105
# 2024

################################################################################

library(tidyverse)
library(readxl)
library(data.table)
library(ggpubr)
library(viridis)

rm(list = ls())
setwd("E:/Research_Projects/MEG_CamCAN/TDE_HMM/output/3_Neurocognitive_analysis")

# DATA WRANGLING ----------------------------------------------------------

age <- rio::import("../cog_data_education.csv")$Age_Cog

# Check if age is positively correlated with each LV. 
# For ease of interpretation all LV will made to positively correlate with age

BSR_temporal <- rio::import("./output_PLS_temporal/BSR_tempo_COG.txt")
BSR_spectral <- rio::import("./output_PLS_spectral/BSR_spectral_COG.txt")
BSR_spectral$V1 <- BSR_spectral$V1*(-1)
BSR_transition <- rio::import("./output_PLS_transition/BSR_trans_COG.txt")
BSR_transition$V1 <- BSR_transition$V1*(-1)


design_names <- c(
  "Stabilizes",
  "Accelerates",
  "Cattell Task",
  "Proverb Task",
  "Naming Task",
  "ToT",
  "Hotel Task",
  "Sentence Comp",
  "Story Recall Task",
  "Verbal Fluency Task",
  "Education"
)


# Temporal
# As mentionned in plot_latent_scores, we do *(-1) to stay consistent
plot_BSR <- cbind(score = BSR_temporal$V1, design_names) %>%
  as.data.frame()
plot_BSR <- plot_BSR %>% mutate_at(vars(score), funs(as.numeric(.)*(-1))) %>% filter(abs(score)>=3)
plot_BSR$design_names <- factor(plot_BSR$design_names) %>%
  fct_reorder(plot_BSR$score, .desc = F)

# remotes::install_github("hrbrmstr/ggchicklet")
library(ggchicklet)
plot_BSR_temporal <- ggplot(
  plot_BSR,
  aes(x = design_names, y = score)
) +
  geom_chicklet(radius = grid::unit(5, "mm"), aes(fill = score, alpha = .8)) +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  scale_y_continuous(limits = c(-20, 21), breaks = seq(-15, 15, 5)) +
  coord_flip() +
  geom_text(aes(y = -0.1, label = design_names, family = 'Arial', fontface = "bold", hjust = -0.1), size = 5) +
  geom_text(aes(y = -0.1, label = score, family = 'Arial', fontface = "bold", hjust = "right"), size = 5) +
  theme_pubr(
    base_size = 16,
    legend = "none",
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_blank()
  )
plot_BSR_temporal


# Spectral
plot_BSR <- cbind(score = BSR_spectral$V1, design_names) %>%
  as.data.frame() 
plot_BSR <- plot_BSR %>% mutate_at(vars(score), funs(as.numeric(.)*(-1))) %>% filter(abs(score)>=3)
plot_BSR$design_names <- factor(plot_BSR$design_names) %>%
  fct_reorder(plot_BSR$score, .desc = F)

plot_BSR_spectral <- ggplot(
  plot_BSR,
  aes(x = design_names, y = score)
) +
  geom_chicklet(radius = grid::unit(5, "mm"), aes(fill = score, alpha = .8)) +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  scale_y_continuous(limits = c(-25, 25), breaks = seq(-25, 25, 5)) +
  coord_flip() +
  geom_text(aes(y = -0.1, label = design_names, family = 'Arial', fontface = "bold", hjust = "right"), size = 5) +
  geom_text(aes(y = -0.1, label = score, family = 'Arial', fontface = "bold", hjust = -.1), size = 5) +
  theme_pubr(
    base_size = 16,
    legend = "none",
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_blank()
  )

plot_BSR_spectral

# Transition
plot_BSR <- cbind(score = BSR_transition$V1, design_names) %>%
  as.data.frame() 
plot_BSR <- plot_BSR %>% mutate_at(vars(score), funs(as.numeric(.)*(-1))) %>% filter(abs(score)>=3)
plot_BSR$design_names <- factor(plot_BSR$design_names) %>%
  fct_reorder(plot_BSR$score, .desc = F)

plot_BSR_transition <- ggplot(
  plot_BSR,
  aes(x = design_names, y = score)
) +
  geom_chicklet(radius = grid::unit(5, "mm"), aes(fill = score, alpha = .8)) +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  scale_y_continuous(limits = c(-25, 25), breaks = seq(-25, 25, 5)) +
  coord_flip() +
  geom_text(aes(y = -0.1, label = design_names, family = 'Arial', fontface = "bold", hjust = "right"), size = 5) +
  geom_text(aes(y = -0.1, label = score, family = 'Arial', fontface = "bold", hjust = -.1), size = 5) +
  theme_pubr(
    base_size = 16,
    legend = "none",
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_blank()
  )

plot_BSR_transition

Rmisc::multiplot(plot_BSR_temporal, plot_BSR_spectral, plot_BSR_transition, cols = 3)


# Average
average <- seq(8)
for (cog_score in seq(8)){
  average[cog_score] <- mean(c(
    abs(BSR_temporal[2+cog_score,1]),
    abs(BSR_spectral[2+cog_score,1]),
    abs(BSR_transition[2+cog_score,1])
  )) %>% as.data.frame()
}

colnames(average) <- c(
  "Cattell Task",
  "Proverb Task",
  "Naming Task",
  "ToT",
  "Hotel Task",
  "Sentence Comp",
  "Story Recall Task",
  "Verbal Fluency Task")

print(average)



# Temporal - FO/LT/INT/SR

BSR_temporal_metrics <- rio::import("E:/Research_Projects/MEG_CamCAN/TDE_HMM/output/3_Neurocognitive_analysis/output_PLS_temporal/BSR_tempo_BRAIN.txt")

metrics_names <- 
  c('fo 1',
    'fo 2',
    'fo 3',
    'fo 4',
    'fo 5',
    'fo 6',
    'fo 7',
    'fo 8',
    'lt 1',
    'lt 2',
    'lt 3',
    'lt 4',
    'lt 5',
    'lt 6',
    'lt 7',
    'lt 8',
    'intv 1',
    'intv 2',
    'intv 3',
    'intv 4',
    'intv 5',
    'intv 6',
    'intv 7',
    'intv 8',
    'sr 1',
    'sr 2',
    'sr 3',
    'sr 4',
    'sr 5',
    'sr 6',
    'sr 7',
    'sr 8'
  )


# These are positively correlated with age
plot_BSR <- cbind(score = BSR_temporal_metrics$V1, metrics_names) %>%
  as.data.frame() 
plot_BSR <- plot_BSR %>% mutate_at(vars(score), funs(as.numeric(.)))
plot_BSR$metrics_names <- factor(plot_BSR$metrics_names) %>%
  fct_reorder(plot_BSR$score, .desc = F)

ggplot(
  plot_BSR,
  aes(x = metrics_names, y = score)
) +
  geom_chicklet(radius = grid::unit(5, "mm"), aes(fill = score, alpha = .8)) +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  scale_y_continuous(limits = c(-20, 15), breaks = seq(-15, 15, 5)) +
  coord_flip() +
  geom_text(aes(y = -0.1, label = metrics_names, family = 'Arial', fontface = "bold", hjust = "right"), size = 5) +
  geom_text(aes(y = -0.1, label = score, family = 'Arial', fontface = "bold", hjust = -1), size = 5) +
  theme_pubr(
    base_size = 16,
    legend = "none",
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_blank()
  )



plot_temporal_metrics <- function(i){
  library(tidyverse)
  library(stringr)
  library(scales)
  library(prismatic)
  
  toplot <- plot_BSR[c(i,i+8,i+8*2,i+8*3),]
  toplot$metrics_names <- factor(toplot$metrics_names, 
                                 levels = unique(toplot$metrics_names))
  levels(toplot$metrics_names) <- c("Total time", "Lifetime", "Interval time", "Switching Rate")
  
  plot_alpha <- ggplot(toplot) +
    
    #make custom panel grid
    geom_hline(yintercept = min(toplot$score), color = "lightgrey") +
    geom_hline(yintercept = quantile(toplot$score)[2], color = "lightgrey") +
    geom_hline(yintercept = 0, color = "gray12", size = 2, linetype = 4) +
    geom_hline(yintercept = quantile(toplot$score)[3], color = "lightgrey") +
    geom_hline(yintercept = max(toplot$score), color = "lightgrey") +
    
    geom_col(aes(
      x = metrics_names, 
      y = score,#is numeric
      fill = score), 
      position = "dodge2",
      show.legend = TRUE,
      alpha = .9) +
    
    if (max(toplot$score) <= 0){ # If only negative
      scale_fill_distiller(palette = "Blues")
    } else if (min(toplot$score) >= 0){ # If only positive
      scale_fill_distiller(palette = "Reds", direction = 1)
    } else if (min(toplot$score) <= 0 && max(toplot$score) >= 0){ # If both
      scale_fill_gradientn("Power (a.u.)",
                           colours = c("#0571B0FF", "#92C5DEFF", "#F4A582FF", "#CA0020FF"))
    }
  
  
  plot_alpha + geom_point(aes(x = metrics_names,
                              y = score),
                          size = 3,
                          color = "gray12") +
    
    #lollipop shaft
    geom_segment(aes(
      x = metrics_names,
      y = min(score),
      xend = metrics_names,
      yend = max(score)),
      linetype = "dashed",
      color = "gray12") +
    
    #theming
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(color = "gray12",
                                     size = 16),
          panel.background = element_rect(fill = "white",
                                          color = "white"),
          panel.grid = element_blank(),
          panel.grid.major.x = element_blank(),
          text = element_text(color = "gray12",
                              family = "Arial"))
  
  # ggsave(path = r"(./results/08_states_networks)",
  #        filename = sprintf("State_%d_RSN_composition.png",i),
  #        width = 20,
  #        height = 20,
  #        units = "cm",
  #        type = "cairo-png",
  #        dpi = 400)


}


# Pattern #1 - DMN
lapply(1, plot_temporal_metrics)
lapply(2, plot_temporal_metrics)
lapply(7, plot_temporal_metrics)
lapply(8, plot_temporal_metrics)

# Pattern #2 - SMN
lapply(4, plot_temporal_metrics)
lapply(5, plot_temporal_metrics)




# LC2

# When plotting the latent scores, we inverted the sign for the connectivity cognitive variable to stay consistent
# So we do the same when examining the BSR

# Temporal
plot_BSR <- cbind(score = BSR_temporal$V2, design_names) %>%
  as.data.frame()
plot_BSR <- plot_BSR %>% mutate_at(vars(score), funs(as.numeric(.))) %>% filter(abs(score)>=3)
plot_BSR$design_names <- factor(plot_BSR$design_names) %>%
  fct_reorder(plot_BSR$score, .desc = F)

# remotes::install_github("hrbrmstr/ggchicklet")
library(ggchicklet)
plot_BSR_temporal <- ggplot(
  plot_BSR,
  aes(x = design_names, y = score)
) +
  geom_chicklet(radius = grid::unit(5, "mm"), aes(fill = score, alpha = .8)) +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  scale_y_continuous(limits = c(-20, 21), breaks = seq(-15, 15, 5)) +
  coord_flip() +
  geom_text(aes(y = -0.1, label = design_names, family = 'Arial', fontface = "bold", hjust = -0.1), size = 5) +
  geom_text(aes(y = -0.1, label = score, family = 'Arial', fontface = "bold", hjust = "right"), size = 5) +
  theme_pubr(
    base_size = 16,
    legend = "none",
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_blank()
  )
plot_BSR_temporal


# Spectral
plot_BSR <- cbind(score = BSR_spectral$V2, design_names) %>%
  as.data.frame() 
plot_BSR <- plot_BSR %>% mutate_at(vars(score), funs(as.numeric(.))) %>% filter(abs(score)>=3)
plot_BSR$design_names <- factor(plot_BSR$design_names) %>%
  fct_reorder(plot_BSR$score, .desc = F)

plot_BSR_spectral <- ggplot(
  plot_BSR,
  aes(x = design_names, y = score)
) +
  geom_chicklet(radius = grid::unit(5, "mm"), aes(fill = score, alpha = .8)) +
  scale_fill_distiller(palette = "Reds", direction = -1) +
  scale_y_continuous(limits = c(-25, 25), breaks = seq(-25, 25, 5)) +
  coord_flip() +
  geom_text(aes(y = -0.1, label = design_names, family = 'Arial', fontface = "bold", hjust = "right"), size = 5) +
  geom_text(aes(y = -0.1, label = score, family = 'Arial', fontface = "bold", hjust = -.1), size = 5) +
  theme_pubr(
    base_size = 16,
    legend = "none",
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_blank()
  )

plot_BSR_spectral


Rmisc::multiplot(plot_BSR_temporal, plot_BSR_spectral, cols = 2)
