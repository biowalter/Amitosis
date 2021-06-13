# Description -------------------------------------------------------------
# Authors: Vitali et al. (2021)
# Figure S1
# Compare empirical IES retention levels between bulk DNA-seq & scDNA-seq
# .........................................................................

# Libraries ---------------------------------------------------------------
library(ggpubr)
library(ggsci)
library(rstatix)
library(tidyverse)

# Import data -------------------------------------------------------------
IES_tab <- read.table("./infiles/amitosis_IES_table.txt", 
                      header = T, stringsAsFactors = F, sep = "\t", quote = "")

# Perform comparison ------------------------------------------------------
# select columns of interest
Test_set_df <- IES_tab %>%
  filter(!is.na(IES_ID) & # rm NAs
           !duplicated(IES_ID)) %>% # drop duplicates
  filter_at(vars(starts_with("IES_COV_P25_sc"), -IES_COV_P25_sc2_D14), 
            all_vars(. > 20)) %>% # filter subset of columns that all fulfill condition
  select(-contains("P25_sc2_D14")) %>%
  select("IES_ID", # select columns of interest
         starts_with("IRS_P25") & ends_with("D7"), 
         starts_with("IRS_P25") & ends_with("D5"),
         starts_with("IRS_P25") & ends_with("D10"),
         starts_with("IRS_P25") & ends_with("D14") ) %>% 
  filter( IRS_P25_sc1_D5 > 0.1 & 
            IRS_P25_sc2_D5 > 0.1 & 
            IRS_P25_sc3_D5 > 0.1 & 
            IRS_P25_sc4_D5 > 0.1)

# reshape track set
Test_set_df_long <- pivot_longer(data = Test_set_df,
                                 cols = matches("IRS"),
                                 values_to = "IRS",
                                 names_to = c("Sample", "Day"),
                                 names_sep = "_",
                                 names_prefix = "IRS_P25_")

Test_set_df_sc_long <- Test_set_df_long %>%
  filter(Sample != "2b") %>% 
  group_by(IES_ID, Day) %>% 
  summarise(IRS = mean(IRS, na.rm = TRUE)) %>% 
  mutate(Sample = "scDNA") %>% 
  relocate(Sample, .after = Day)

Test_set_df_mc_long <- Test_set_df_long %>% 
  filter(Sample == "2b")

Test_set_df_combo_long <- bind_rows(Test_set_df_mc_long, Test_set_df_sc_long)
Test_set_df_combo_long$Day <- factor(Test_set_df_combo_long$Day, 
                                     levels = c("D7", "D5", "D10", "D14"))

# rm IRS = 1.0
Test_set_df_combo_long <- Test_set_df_combo_long %>% 
  filter(IRS != 1.0)

# Rename 2b to mcDNA
Test_set_df_combo_long$Sample[Test_set_df_combo_long$Sample=="2b"] <- "mcDNA"

# compute summary statistics
Test_set_df_combo_long %>% 
  group_by(Day) %>%
  get_summary_stats(IRS, type = "median_iqr")

# computation
stat.test <- Test_set_df_combo_long %>% 
  wilcox_test(IRS ~ Day, alternative = "two.sided", 
              paired = FALSE, ref.group = "D7",
              p.adjust.method = "BH") %>%
  add_significance()
stat.test

# add significance to plot
# compare IRS (e.g. mcDNA vs scDNA)
ggboxplot(Test_set_df_combo_long, x= "Day", y = "IRS",
          add = c("jitter"),
          bxp.errorbar = T,
          notch = T, color = "Sample",
          alpha = 0.1,
          palette = c("dark orange", "#00468BFF"),
          width = 0.5) +
  labs(title="", x="Day", y = "Retention levels (IRS)") +
  theme_bw(base_size = 25) +
  theme(legend.title = element_text(color = "black", size = 20, face = "plain"),
        axis.text = element_text(color = "black"),
        axis.title.x = element_text(color = "black", size = 25, face = "plain", vjust=-1),
        axis.title.y = element_text(color = "black", size = 25, face = "plain", vjust= 1, 
                                    angle = 90, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.position = "top", plot.margin = unit(c(0,4,1,2), "cm")) + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif",
                     y.position = c(1.10, 1.20, 1.30)) +
  scale_y_continuous(breaks = seq(0.0, 1.0, by = 0.25))

# sample size
Test_set_df_combo_long %>% 
  group_by(Day) %>% 
  summarise(n())
