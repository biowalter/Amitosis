# Description -------------------------------------------------------------
# Authors: Vitali et al. (2021)
# Random measurement error of IES Retention Levels
# .........................................................................

# Libraries ---------------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(ggsci)
library(RColorBrewer)
library(patchwork)

# Load data ---------------------------------------------------------------
## Full IES Track Set data frame
Track_set_df_red<- read_csv("./infiles/track_set_df_red_chromo_haplo_no_dist")
spec(Track_set_df_red)

# Measurement Error of IES Retention Levels --------------------------------
## Empirical error on IRS measurements

## load MIRET output with read.tab_rec
load("./fun/read.table_recursive_v2.RData")
IRS_error <- read.table_recursive(mypath = "./infiles/MIRET/", mypattern = "MIRET", sep = "\t", h = T, stringsAsFactors = F)
sample_names <- lapply(names(IRS_error), function(x) paste(unlist(strsplit(x, split = "_", fixed = T))[2:3], collapse = "_"))
sample_names <- unlist(lapply(sample_names, function(x) unlist(strsplit(x, split = ".", fixed = T))[1]))

## Add names to data frames
for (i in 1:length(sample_names)) {
  IRS_error[[i]]$sample <- sample_names[i]
}

IRS_error <- do.call(rbind, IRS_error)

## filter out unsupported IESs (set new filter for boundary method)
IRS_error <- IRS_error[(IRS_error$SUPPORT_MAC + IRS_error$SUPPORT_LEFT + IRS_error$SUPPORT_RIGHT > 20), ]

## add support column
IRS_error$SUPPORT_MAC_RIGHT <- (IRS_error$SUPPORT_MAC + IRS_error$SUPPORT_RIGHT)

## rm unsupported IESs
IRS_error <- IRS_error[IRS_error$SUPPORT_MAC_RIGHT > 20, ]

# Test difference btw left and right IRS (binomial test) ------------------
## ctr score = left IRS
## successes = SUPPORT_RIGHT
## trials = SUPPORT_MAC + SUPPORT_RIGHT

## define binomial test function
binom_test_tab <- function(tab){
  if(is.na(tab[1]) | (tab[2]==0 & tab[3]==0)){return (NA)}
  else{return (binom.test(x=tab[2], n=tab[3], p=tab[1], alternative="two.sided")$p.value)}
}

## Apply function to binomial table
pvalues <- apply(IRS_error[, c(13, 18, 22)], 1, binom_test_tab)
IRS_error$padj <- p.adjust(pvalues, method="BH")

## discard IESs with significant difference between lx & rx boundary
### filter somatic IESs
### count before and after removing putative events of differential boundary usage
as_tibble(IRS_error) %>% 
  filter(!is.na(ID), 
         RETENTION_SCORE_LEFT != 1.0, 
         RETENTION_SCORE_LEFT >= 0.1, 
         RETENTION_SCORE_RIGHT >= 0.1,
         padj > 0.05) %>% # toggle
  summarise(n())

IRS_error <- IRS_error[IRS_error$padj > 0.05 & 
                         IRS_error$RETENTION_SCORE_LEFT != 1.0 &
                         (IRS_error$RETENTION_SCORE_LEFT >= 0.1 & IRS_error$RETENTION_SCORE_RIGHT >= 0.1), ]
### rm IESs shorter than 50 bp
#IRS_error <- IRS_error[nchar(IRS_error$SEQUENCE) > 50, ]
#View(IRS_error) # inspect df
### compute error (cv) and mean boundary scores
IRS_error$SD <- apply(select(IRS_error, matches("retention")), 1, sd)
IRS_error$IRS_mean <-  apply(select(IRS_error, matches("retention")), 1, mean)
IRS_error <- IRS_error[IRS_error$SD > 0, ] # rm sd = 0

# Plot - IRS_SD / IRS_mean ~ IRS_mean -----------------------------------
relative_error <- ggplot() + 
  geom_point(mapping = aes(x = IRS_mean, y = SD / IRS_mean), 
             data = IRS_error, 
             alpha = 0.2, size = 2, color="dark blue") +
  geom_smooth(aes(x = IRS_mean, y = SD / IRS_mean), 
              data = IRS_error, color="dark red",
              method = "loess", size = 1) +
  theme_bw(base_size = 20)

# add bins category to df
IRS_error <- IRS_error %>%
  mutate(bin = cut_width(IRS_mean,  width = 0.05, boundary = 0))


# Plot - Random Error vs observed variation in IRSs -----------------------
## violin plot

## combine IRS_error and Track_set_df_red
### subset and add label
error <- as_tibble(IRS_error) %>% 
  select(SD) %>%
  mutate(set = "error")
### compute sd across REPs for D14
Track_set_df_red$sd_obs <- Track_set_df_red %>%
  select(ends_with("D14")) %>%
  apply(1, sd)
### subset and add label
observed <- Track_set_df_red %>%
  select(sd_obs) %>%
  rename(SD = sd_obs) %>%
  mutate(set = "observed")
### bind dfs
combo_sd <- bind_rows(error, observed)

## violin plot
ggplot(combo_sd, aes(x = set, y = SD)) +
  geom_violin(aes(fill = set), trim = TRUE, color = "black") +
  scale_fill_brewer(palette="Blues") +
  geom_boxplot(width=0.1, fill = "white") + 
  theme_bw(base_size = 20)

## Perform pairwise comparisons
compare_means(SD ~ set,  data = combo_sd)

## compute summary statistics
combo_sd %>% 
  group_by(set) %>%
  get_summary_stats(SD, type = "median_iqr")

## computation
stat.test <- combo_sd %>%
  wilcox_test(SD ~ set, alternative = "less", paired = FALSE) %>%
  add_significance()
stat.test

## effect size
combo_sd %>% 
  wilcox_effsize(SD ~ set, alternative = "less", paired = FALSE)

## violin plot with stats (Fig5D)
violin <- ggpubr::ggviolin(combo_sd, x = "set", y = "SD", 
                           trim = TRUE, fill = "set", palette = c("#D8988B", "#91D1C2")) +
  geom_boxplot(width=0.1, fill = "white") +
  coord_cartesian(ylim=c(0, 0.25)) +
  stat_compare_means(comparisons = list(c("error", "observed"))) +
  #stat_compare_means(label.y = 0.28) +
  xlab("IES set") +
  theme_bw(base_size = 30) +
  theme(legend.position="top",
        panel.background = element_rect(fill = 'white smoke')) +
  rremove("legend")
violin

# Plot - IRS_SD / IRS_mean ~ IRS_mean (binned) ----------------------------
## 0.05-wide IRS_mean bins
pal <- colorRampPalette(pal_npg("nrc", alpha = 0.5)(7))(18)

relative_error_binned <- ggplot(data = IRS_error, mapping = aes(x = IRS_mean, y = SD / IRS_mean)) +
  geom_boxplot(aes(group = bin, fill = bin), show.legend = FALSE) + 
  scale_fill_manual(values = pal) +
  xlab("Mean boundary IRS (bIRS, binned)") +
  ylab(expression('SD / bIRS (CV)')) +
  theme_bw(base_size = 20) +
  theme(legend.position="top",
        panel.background = element_rect(fill = 'white smoke'),
        plot.margin = margin(t = 1.5, r = 1, b = 1, l = 0, unit = "cm"))
relative_error_binned

## random error estimates for IRS measurements (binned, width = 0.05, n = 18)
IRS_error %>% 
  group_by(cut_width(IRS_mean, width = 0.05, boundary = 0)) %>%
  summarise(median(SD))

## distribution of absolute and relative measurement errors for IRSs
summary(IRS_error$SD) # absolute
summary(IRS_error$SD / IRS_error$IRS_mean) # relative


# Plot - IRS_SD / IRS_mean ~ IRS_mean Error vs. Observed ------------------
relative_error1 <- ggplot() + 
  geom_point(mapping = aes(x = IRS_mean, y = SD / IRS_mean, color="dark gray"), 
             data = IRS_error, # error
             alpha = 0.5, size = 2) +
  geom_smooth(aes(x = IRS_mean, y = SD / IRS_mean, color="dark red"), 
              data = IRS_error, # fit to error
              method = "loess", fill = "light salmon") +
  geom_point(mapping = aes(x=IRS_P25_sc_D5, y=sd_obs / IRS_P25_sc_D5, color="blue"), 
             data = Track_set_df_red, # empirical
             alpha = 0.2, size = 2) + 
  geom_smooth(aes(x=IRS_P25_sc_D5, y=sd_obs / IRS_P25_sc_D5, color="dark blue"), 
              data = Track_set_df_red, # fit to empirical
              method = "loess", fill = "light salmon") +
  scale_color_identity(name = "set",
                       breaks = c("dark red", "dark blue"),
                       labels = c("fit error", "fit observed"),
                       guide = "legend") +
  labs(x = expression('Mean IRS (IRS'['0']*', gen=17)'),
       y = expression('SD / IRS'['0']*' (CV, gen=49)')) +
  theme_bw(base_size = 20) +
  theme(legend.position="top",
        panel.background = element_rect(fill = 'white smoke'),
        plot.margin = margin(t = 0, r = 0, b = 1, l = 1, unit = "cm"))
relative_error1

# Arrange error plots
ggarrange(relative_error_binned, relative_error1, violin, nrow = 3, 
          align = "v", labels = c("A", "B", "C"), label.y = c(1, 1.1, 1))

relative_error_binned + relative_error1 + plot_annotation(tag_levels = 'A')

# Plot - Change of SD / IRS_mean with IRS0 --------------------------------
## Expected vs observed (empirical SD / IRS_mean at D14)
## compute sd across REPs for D14
Fig5C <- ggplot() +
  geom_point(mapping = aes(x = IRS_mean, y = SD / IRS_mean, color="dark gray"), 
             data = IRS_error, # rel error
             alpha = 0.2, size = 2) +
  geom_smooth(aes(x = IRS_mean, y = SD / IRS_mean), 
              data = IRS_error, color="dark red", # fit to rel error
              method = "loess", fill = "light salmon") +
  geom_point(mapping = aes(x=IRS_P25_sc_D5, y=sd_obs / IRS_P25_sc_D5, color="blue"), 
             data = Track_set_df_red, # empirical
             alpha = 0.4, size = 2) + 
  geom_smooth(aes(x=IRS_P25_sc_D5, y=sd_obs / IRS_P25_sc_D5), 
              data = Track_set_df_red, color="dark blue", # fit to empirical
              method = "loess") +
  geom_point(mapping = aes(x=IRS_P25_sc_D5, y=sd_chromo / IRS_P25_sc_D5, color="dark orange"), 
             data = Track_set_df_red,  # simulated chromosomal
             alpha = 0.8, size = 2) +
  geom_point(mapping = aes(x=IRS_P25_sc_D5, y=sd_haplo / IRS_P25_sc_D5, color="#079787"), 
             data = Track_set_df_red, # simulated haploid
             alpha = 0.8, size = 2) +
  scale_color_identity(name = "",
                       breaks = c("dark gray", "blue", "dark orange", "#079787"),
                       labels = c("error", "observed", "predicted chrom.", "predicted haplo."),
                       guide = "legend") +
  labs(x = expression('Mean IRS (IRS'['0']*', gen=17)'),
       y = expression('SD / IRS'['0']*' (CV, gen=49)')) +
  theme_bw(base_size = 25) +
  theme(plot.background=element_rect(fill = "white"),
        axis.title.x = element_text(color = "black", size = 25, face = "plain", vjust=-2),
        axis.title.y = element_text(color = "black", size = 25, face = "plain", vjust= 1.0, angle = 90, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        panel.background = element_rect(fill = 'white smoke'),
        legend.position="top", legend.box = "horizontal", 
        legend.spacing.x = unit(0.1, 'cm'),
        plot.margin = unit(c(1,3,1,1), "cm"))

Fig5C <- ggpar(Fig5C, 
               legend = "top",
               font.legend = c(20, "plain", "black"))

## Run script "SA_obs_vs_exp.R" to combine figures
#ggarrange(CI_haplo, CI_chromo, Fig5C, labels = "AUTO", align = "v",
          #font.label = list(size = 30, face = "plain", color ="black"))
