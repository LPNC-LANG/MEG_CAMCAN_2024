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
library(flexplot)

rm(list = ls())
setwd("E:/Research_Projects/MEG_CamCAN/TDE_HMM/output/3_Neurocognitive_analysis")

# DATA WRANGLING ----------------------------------------------------------

age <- rio::import("./cog_data_education.csv")$Age_Cog

# Check if age is positively correlated with each LV. 
# For ease of interpretation all LV will made to positively correlate with age
# library(flexplot)
# e.g., flexplot::flexplot(VD~age, data = cbind(VD = Lx_temporal$V1,age) %>% as.data.frame())
Lx_temporal <- rio::import("./output_PLS_temporal/Lx.csv")
Ly_temporal <- rio::import("./output_PLS_temporal/Ly.csv")

Lx_spectral <- rio::import("./output_PLS_spectral/Lx.csv")
Lx_spectral$V1 <- Lx_spectral$V1*(-1)
Ly_spectral <- rio::import("./output_PLS_spectral/Ly.csv")
Ly_spectral$V1 <- Ly_spectral$V1*(-1)

Lx_transition <- rio::import("./output_PLS_transition/Lx.csv")

Lx_transition$V1 <- Lx_transition$V1*(-1)
Ly_transition <- rio::import("./output_PLS_transition/Ly.csv")
Ly_transition$V1 <- Ly_transition$V1*(-1)

# DATA ANALYSIS -----------------------------------------------------------

model_temporal <- cbind(
  Lx1 = Lx_temporal[,1],
  Lx2 = Lx_temporal[,2],
  Ly1 = Ly_temporal[,1],
  Ly2 = Ly_temporal[,2]
) %>% as.data.frame()

ggplot_temporal <- ggplot(model_temporal, aes(x = Lx1, y = Ly1, color = age)) +
  # geom_smooth(method = "gam", color = "black", size = 2, alpha = 0.3, formula = y ~ s(x, k=3)) +
  geom_point(aes(color = age), size = 5) +
  # scale_y_continuous(limits = c(-2, 4), breaks = seq(-2, 4, 2)) +
  scale_color_viridis(name = 'Cognitive age', discrete = F, option = "magma")+
  labs(
    color = "Age groups",
    title = "LC1 - Temporal"
  ) +
  theme_pubr(base_size = 14)

ggplot_temporal

latent_scores <- cbind(model_temporal_brain = model_temporal$Lx1, 
                       model_temporal_cog = model_temporal$Ly1, age = age) %>% as.data.frame() %>% 
  pivot_longer(c("model_temporal_brain", "model_temporal_cog"), names_to = "PLS model", values_to = "Values")

latent_scores %>%
  group_by(`PLS model`) %>% 
  mutate(Values = scale(Values)) %>% 
  ggplot(aes(age, Values, color = `PLS model`)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_jitter(alpha=0.1) +
  geom_smooth(linewidth = 3, method = "loess",alpha = 0.2) +
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  # scale_y_continuous(breaks = seq(-1.5, 1, 0.5)) +
  coord_cartesian(ylim = c(-3, 3)) +
  scale_color_manual(values = c("#0571B0FF", "#CA0020FF")) +
  theme_pubclean(base_size = 16) +
  labs(
    y = "Latent score",
    x = "Age (years)"
  )



model_spectral <- cbind(
  Lx1 = Lx_spectral[,1],
  Lx2 = Lx_spectral[,2],
  Ly1 = Ly_spectral[,1],
  Ly2 = Ly_spectral[,2]
) %>% scale() %>% as.data.frame()

ggplot_spectral <- ggplot(model_spectral, aes(x = Lx1, y = Ly1, color = age)) +
  geom_point(aes(color = age), size = 5) +
  scale_color_viridis(name = 'Cognitive age', discrete = F, option = "magma")+
  theme_pubr(base_size = 14)

ggplot_spectral

latent_scores <- cbind(model_spectral_brain = model_spectral$Lx1, 
                       model_spectral_cog = model_spectral$Ly1, age = age) %>% as.data.frame() %>% 
  pivot_longer(c("model_spectral_brain", "model_spectral_cog"), names_to = "PLS model", values_to = "Values")

latent_scores %>%
  group_by(`PLS model`) %>% 
  mutate(Values = scale(Values)) %>% 
  ggplot(aes(age, Values, color = `PLS model`)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_jitter(alpha=0.1) +
  geom_smooth(linewidth = 3, method = "loess",alpha = 0.2) +
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  coord_cartesian(ylim = c(-2, 2)) +
  scale_color_manual(values = c("#0571B0FF", "#CA0020FF")) +
  theme_pubclean(base_size = 16) +
  labs(
    y = "Latent scores",
    x = "Age (years)"
  )


model_transition <- cbind(
  Lx1 = Lx_transition[,1],
  Lx2 = Lx_transition[,2],
  Ly1 = Ly_transition[,1],
  Ly2 = Ly_transition[,2]
) %>% as.data.frame()

ggplot_transition <- ggplot(model_transition, aes(x = Lx1, y = Ly1, color = age)) +
  # geom_smooth(method = "gam", color = "black", size = 2, alpha = 0.3, formula = y ~ s(x, k=3)) +
  geom_point(aes(color = age), size = 5) +
  # scale_y_continuous(limits = c(-2, 4), breaks = seq(-2, 4, 2)) +
  scale_color_viridis(name = 'Cognitive age', discrete = F, option = "magma")+
  labs(
    color = "Age groups",
    title = "LC1 - transition"
  ) +
  theme_pubr(base_size = 14)

ggplot_transition

latent_scores <- cbind(model_transition_brain = model_transition$Lx1, 
                       model_transition_cog = model_transition$Ly1, age = age) %>% as.data.frame() %>% 
  pivot_longer(c("model_transition_brain", "model_transition_cog"), names_to = "PLS model", values_to = "Values")

latent_scores %>%
  group_by(`PLS model`) %>% 
  mutate(Values = scale(Values)) %>% 
  ggplot(aes(age, Values, color = `PLS model`)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_jitter(alpha=0.1) +
  geom_smooth(linewidth = 3, method = "loess",alpha = 0.2) +
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  # scale_y_continuous(breaks = seq(-1.5, 1, 0.5)) +
  coord_cartesian(ylim = c(-3, 3)) +
  scale_color_manual(values = c("#0571B0FF", "#CA0020FF")) +
  theme_pubclean(base_size = 16) +
  labs(
    y = "Latent score",
    x = "Age (years)"
  )

Rmisc::multiplot(ggplot_temporal, ggplot_spectral, ggplot_transition, cols=3)


# Combining both latent components ----
library(ggpubr)

# We want to plot a decrease across the lifespan
# When interpreting BSR values, we will need to do *(-1) as well

latent_scores_cog <- cbind(Ly1_tempo = model_temporal$Ly1*(-1), 
                           # Ly1_spectral = model_spectral$Ly1*(-1),
                           Ly1_transition = model_transition$Ly1*(-1),
                           age = age) %>% as.data.frame() %>% 
  pivot_longer(c("Ly1_tempo", 
                 # "Ly1_spectral", 
                 "Ly1_transition"
  ), 
  names_to = "Domain of analysis", values_to = "Values")

ggplot_cog <- latent_scores_cog %>%
  group_by(`Domain of analysis`) %>% 
  mutate(Values = scale(Values)) %>% 
  ggplot(aes(age, Values, color = `Domain of analysis`)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_jitter(alpha=0.1, size = 3) +
  geom_smooth(linewidth = 2, method = "gam", alpha = 0.2, formula = y~s(x, k = 3)) +
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  coord_cartesian(ylim = c(-2,2)) +
  scale_color_manual(values = c(
    "#0571B0FF", "darkblue"
    # "#CA0020FF",
  ),
  labels = c(
    # "Spectral", 
    "Temporal", 
    "Transition")) +
  theme_pubr(base_size = 14,
             legend = "none"
  ) +
  labs(
    y = "Latent cognitive scores (Component 1)",
    x = "Age (years)"
  ) 
# Temporal
# annotate("text", 
#          x = 85, y = -0.8,
#          label = "edf = 1.91") +
# # Spectral
# annotate("text", 
#          x = 85, y = -1.1,
#          label = "edf = 1.96")



mgcv::gam(Ly1 ~ s(age, k = 3), data = model_temporal) %>% summary()
mgcv::gam(Ly1 ~ s(age, k = 3), data = model_spectral) %>% summary()
mgcv::gam(Ly1 ~ s(age, k = 3), data = model_transition) %>% summary()


# We want to plot an across the lifespan
latent_scores_brain <- cbind(Lx1_tempo = model_temporal$Lx1, 
                             # Lx1_spectral = model_spectral$Lx1,
                             Lx1_transition = model_transition$Lx1,
                             age = age) %>% as.data.frame() %>% 
  pivot_longer(c("Lx1_tempo", 
                 # "Lx1_spectral", 
                 "Lx1_transition" 
  ),
  names_to = "Domain of analysis", values_to = "Values")

ggplot_brain <- latent_scores_brain %>%
  group_by(`Domain of analysis`) %>% 
  mutate(Values = scale(Values)) %>% 
  ggplot(aes(age, Values, color = `Domain of analysis`)) + # We interpret the latent brain variable in the original sign
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_jitter(alpha=0.1, size = 3) +
  geom_smooth(linewidth = 2, method = "gam", alpha = 0.2, formula = y ~s(x,k=3)) +
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  coord_cartesian(ylim = c(-2,2)) +
  scale_color_manual(values = c(
    "#0571B0FF", "darkblue" 
    # "#CA0020FF"
  ),
  labels = c(
    # "Connectivity", 
    # "Spectral", 
    "Duration & Frequency\n (within states)", "\nState-to-state transitions\n (between states)")) +
  theme_pubr(base_size = 14) +
  theme(legend.position = c(0.75,0.18
  )) +
  labs(
    y = "Latent brain scores (Component 1)",
    x = "Age (years)"
  )

mgcv::gam(Lx1 ~ s(age, k = 3), data = model_temporal) %>% summary()
mgcv::gam(Lx1 ~ s(age, k = 3), data = model_spectral) %>% summary()
mgcv::gam(Lx1 ~ s(age, k = 3), data = model_transition) %>% summary()

Rmisc::multiplot(ggplot_cog, ggplot_brain, cols = 2)



# LC2

latent_scores_cog <- cbind(
  # Ly2_tempo = model_temporal$Ly2,
  Ly2_spectral = model_spectral$Ly2,
  # Ly2_transition = model_transition$Ly2,
  # Ly2_conn = model_conn$Ly2, 
  age = age) %>% as.data.frame() %>% 
  pivot_longer(c(
    # "Ly2_tempo", 
    "Ly2_spectral",
    # "Ly2_transition"
    # "Ly2_conn"
  ), 
  names_to = "Domain of analysis", values_to = "Values")

ggplot_cog <- latent_scores_cog %>%
  group_by(`Domain of analysis`) %>% 
  mutate(Values = scale(Values)) %>% 
  ggplot(aes(age, Values, color = `Domain of analysis`)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_jitter(alpha=0.1, size = 3) +
  geom_smooth(linewidth = 2, method = "gam", alpha = 0.2, formula = y~s(x, k = 3)) +
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  coord_cartesian(ylim = c(-2,2)) +
  scale_color_manual(values = c(
    # "#0571B0FF", "darkblue", 
    "#CA0020FF"
    # "#CA0020FF", "orange"
  ),
  labels = c(
    "Spectral" 
    # "Temporal (within states)", "Temporal (between states)"
  )) +
  theme_pubr(base_size = 14, legend = "none") +
  labs(
    y = "Latent cognitive scores (Component 2)",
    x = "Age (years)"
  ) 
# Temporal
# annotate("text", 
#          x = 85, y = -0.8,
#          label = "edf = 1.91") +
# # Spectral
# annotate("text", 
#          x = 85, y = -1.1,
#          label = "edf = 1.96")

ggplot_cog


mgcv::gam(Ly2 ~ s(age, k = 3), data = model_temporal) %>% summary()
mgcv::gam(Ly2 ~ s(age, k = 3), data = model_spectral) %>% summary()
mgcv::gam(Ly2 ~ s(age, k = 3), data = model_transition) %>% summary()

latent_scores_brain <- cbind(
  # Lx2_tempo = model_temporal$Lx2, 
  Lx2_spectral = model_spectral$Lx2,
  # Lx2_transition = model_transition$Lx2,
  # Lx2_conn = model_conn$Lx2, 
  age = age) %>% as.data.frame() %>% 
  pivot_longer(c(
    # "Lx2_tempo", 
    "Lx2_spectral", 
    # "Lx2_transition" 
    # "Lx2_conn"
  ),
  names_to = "Domain of analysis", values_to = "Values")

ggplot_brain <- latent_scores_brain %>%
  group_by(`Domain of analysis`) %>% 
  mutate(Values = scale(Values)) %>% 
  ggplot(aes(age, Values, color = `Domain of analysis`)) + # We interpret the latent brain variable in the original sign
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_jitter(alpha=0.1, size = 3) +
  geom_smooth(linewidth = 2, method = "gam", alpha = 0.2, formula = y ~s(x,k=3)) +
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  coord_cartesian(ylim = c(-2,2)) +
  scale_color_manual(values = c(
    # "#0571B0FF", "darkblue", 
    "#CA0020FF"
    # , "orange"
  ),
  labels = c(
    "Spectral"
    # "Temporal (within states)", "Temporal (between states)"
  )) +
  theme_pubr(base_size = 14) +
  theme(legend.position = c(0.75,0.09
  )) +
  labs(
    y = "Latent brain scores (Component 2)",
    x = "Age (years)"
  )

mgcv::gam(Lx2 ~ s(age, k = 3), data = model_temporal) %>% summary()
mgcv::gam(Lx2 ~ s(age, k = 3), data = model_spectral) %>% summary()
mgcv::gam(Lx2 ~ s(age, k = 3), data = model_transition) %>% summary()

Rmisc::multiplot(ggplot_cog, ggplot_brain, cols = 2)

