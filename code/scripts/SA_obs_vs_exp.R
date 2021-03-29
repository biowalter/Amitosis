# Description -------------------------------------------------------------
# Authors: Vitali et al. (2021)
# Compare observed and theoretical variation of IES retention levels
# ~50 amitotic divisions post self-fertilization
# .........................................................................

# Libraries ---------------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(rstatix)
library(fBasics)
library(intervcomp)
library(tidyverse)
library(ggsci)
library(RColorBrewer)
library(patchwork)

# Define functions --------------------------------------------------------
## SENES core module
## 95% Confidence Interval
senes_core <- function(k, Chr, N, IRS0, GEN, hetero_df, hetero_df_all) {
  for (w in 1:length(IRS0)) {
    m = round(2*IRS0[w]*k)-1 # number of success states, IES+ copies available (integer)
    input = m / 2 / k # input_ratio
    if (input > 0.5) {
      m = round(2*(1-IRS0[w])*k)+1
    }
    n = N - m # number of failure states, IES- copies available (or non IES+ copies)
    x = seq(0, k*Chr, 1) # All possible successes
    y = dhyper(x, m[1], n, k*Chr) # 1st GEN
    
    for (j in 1:length(GEN)) {
      d <- list()
      for (i in 1:length(which(ifelse((y > 0), TRUE, FALSE)))) {
        d[[i]] <- y[which(ifelse((y > 0), TRUE, FALSE))][i] * 
          dhyper(x, 2*(which(ifelse((y > 0), TRUE, FALSE))[i]-1), N-2*(which(ifelse((y > 0), TRUE, FALSE))[i]-1), k*Chr)
      }
      y <- d[[1]]
      for (i in 2:length(d)) {
        y <- y + d[[i]] 
      }
      x1 = seq(0, k)
      y1 = c(y[1:(k)], sum(y[(k+1):length(y)])) / 2
      # ......................Magic Block
      if (input == 0.5) {
        y2 = y1 + rev(y1)
      }
      else if (input < 0.5) {
        y2 = y1*2
      }
      else if (input > 0.5) {
        y2 = rev(y1)*2
      }
      # .....................End of Magic
      spread <-  sqrt(sum(x1^2 * y2) - sum(x1 * y2)^2) # calculate sd for num of IES+ copies
      hetero_df[j, ] <- cbind(j+1, IRS0[w], y2[1]+y2[k+1], sum(y2[2:(k)]), spread, list(list(y2))) # SD of each iteration
    }
    
    hetero_df_all <- rbind(hetero_df_all, hetero_df)
  }
  hetero_df_all$sd <- as.numeric(hetero_df_all$sd)/k
  return(hetero_df_all)
}
LB <- function(x, SD) {
  if (x < 1) {
    LB = x - 2*SD
    ifelse(LB < 0, 0, LB)
  } else  {1}
}
UB <- function(x, SD) {
  if (x > 0) {
    UB = x + 2*SD
    ifelse(UB > 1, 1, UB)
  }else {0}
}

# Import data -------------------------------------------------------------
IES_tab <- read_tsv("./infiles/amitosis_IES_table.txt")
spec(IES_tab)

# Manipulate data ---------------------------------------------------------
## subset IES table and filter (aggressive)
## IES_COV > 20 in all 11 sc samples
## somatic IESs at D5 (IRS > 0.1)

## select columns of interest
IES_tab1 <- IES_tab %>%
  filter(!is.na(IES_ID) & # rm NAs
           !duplicated(IES_ID)) %>% # drop duplicates
  filter_at(vars(starts_with("IES_COV_P25_sc"), -IES_COV_P25_sc2_D14), 
            all_vars(. > 20)) %>% # filter subset of columns that all fulfill condition
  select(-contains("P25_sc2_D14")) %>%
select("IES_ID", # select columns of interest
         starts_with("IRS_P25") & ends_with("D5"), 
         starts_with("IRS_P25") & ends_with("D10"),
         starts_with("IRS_P25") & ends_with("D14") ) %>%
  filter( IRS_P25_sc1_D5 > 0.1 & 
            IRS_P25_sc2_D5 > 0.1 & 
            IRS_P25_sc3_D5 > 0.1 & 
            IRS_P25_sc4_D5 > 0.1) %>%
  print(n = Inf)

# Compare IRS sd values across time points --------------------------------
## plot sd values
### calculate sd observed
sd_df <- IES_tab1 %>%
  select("IES_ID", ends_with("D5"), ends_with("D10"), ends_with("D14")) %>% 
  mutate(SD_D5 = rowSds(across(ends_with("D5")))) %>% #  transmute(SD_D5 = apply(across(ends_with("D5")), 1 , sd))
  mutate(SD_D10 = rowSds(across(ends_with("D10")))) %>%
  mutate(SD_D14 = rowSds(across(ends_with("D14")))) %>%
  select(IES_ID, matches("SD"))

IES_tab1 %>%
  select("IES_ID", ends_with("D5"), ends_with("D10"), ends_with("D14")) %>% 
  mutate(SD_D5 = across(ends_with("D5"), sd)) %>%
  mutate(SD_D10 = across(ends_with("D10"), sd)) %>%
  mutate(SD_D14 = across(ends_with("D14"), sd)) %>%
  select(IES_ID, matches("SD"))

### reshape to long
sd_df_long <- pivot_longer(data = sd_df,
                           cols = matches("SD"),
                           values_to = "sd_obs",
                           names_to = "day",
                           names_prefix = "SD_")

### summarize median sd by day
sd_df_long %>%
  group_by(day) %>%
  summarise(sd_median = median(sd_obs, na.rm = TRUE),
            sd_mean = mean(sd_obs, na.rm = TRUE))

### compute summary statistics
sd_df_long %>%
  group_by(day) %>%
  get_summary_stats(sd_obs, type = "median_iqr")

### Assumptions and preliminary tests
### differences between paired samples should be distributed symmetrically around the median
sd_df %>%
  mutate(diff = SD_D14 - SD_D5) %>%
  gghistogram(x = "diff", y = "..density..", 
              fill = "steelblue",bins = 10, add_density = TRUE)

### computation
stat.test <- sd_df_long %>%
  wilcox_test(sd_obs ~ day, alternative = "greater", paired = TRUE) %>%
  add_significance()
stat.test

### effect size
sd_df_long  %>%
  wilcox_effsize(sd_obs ~ day, alternative = "greater", paired = TRUE)

### violin plot
sd_df_long$day <- factor(as.factor(sd_df_long$day), levels = c("D5", "D10", "D14"))
violin1 <- ggplot(sd_df_long, aes(x = day, y = sd_obs)) +
  geom_violin(aes(fill = day), trim = TRUE, color = "black") +
  coord_cartesian(ylim=c(0, 0.30)) +
  scale_fill_brewer(palette="Blues") +
  geom_boxplot(width=0.1, fill = "white") + 
  stat_compare_means(comparisons = list(c("D10", "D5"),
                                        c("D14", "D5"), 
                                        c("D14", "D10")),
                     method.args = list(alternative = "greater"), paired = TRUE) +
  theme_bw(base_size = 20) +
  rremove("legend")

## explore change in IRS standard deviation over time using sd ratios 
sd_ratio_df <- IES_tab1 %>%
  select("IES_ID", ends_with("D5"), ends_with("D10"), ends_with("D14")) %>% 
  mutate(SD_D5 = rowSds(across(ends_with("D5")))) %>%
  mutate(SD_D10 = rowSds(across(ends_with("D10")))) %>%
  mutate(SD_D14 = rowSds(across(ends_with("D14")))) %>% 
  mutate(R_D10_to_D5 = SD_D10 / SD_D5) %>%
  mutate(R_D14_to_D5 = SD_D14 / SD_D5) %>%
  mutate(R_D14_to_D10 = SD_D14 / SD_D10) %>%
  select(IES_ID, matches("R_")) %>%
  filter(R_D14_to_D5 != "NaN", R_D14_to_D5 != "Inf") %>%
  print(n = Inf)

sd_ratio_df %>%
  filter(R_D14_to_D5 > 1) %>%
  print(n = Inf)

# Plot SD_ratios ----------------------------------------------------------
## reshape to long
sd_ratio_df_long <- pivot_longer(data = sd_ratio_df,
                                 cols = matches("R"),
                                 values_to = "sd_ratio",
                                 names_to = "time",
                                 names_prefix = "R_")

## compute summary statistics
sd_ratio_df_long %>% 
  group_by(time) %>%
  get_summary_stats(sd_ratio, type = "median_iqr")

## Assumptions and preliminary tests
## differences between paired samples should be distributed symmetrically around the median
sd_ratio_df %>% 
  mutate(diff = R_D14_to_D5 - R_D10_to_D5) %>%
  gghistogram(x = "diff", y = "..density..", 
              fill = "steelblue",bins = 10, add_density = TRUE)

## computation
stat.test1 <- sd_ratio_df_long %>% 
  wilcox_test(sd_ratio ~ time, alternative = "less", paired = TRUE) %>%
  add_significance()
stat.test1

## effect size
sd_ratio_df_long  %>% 
  wilcox_effsize(sd_ratio ~ time, alternative = "less", paired = TRUE)

## violin plot
violin2 <- ggplot(sd_ratio_df_long, aes(x = time, y = sd_ratio)) +
  geom_violin(aes(fill = time), trim = TRUE, color = "black") +
  coord_cartesian(ylim=c(0, 13)) +
  scale_fill_brewer(palette="Blues") +
  geom_boxplot(width=0.1, fill = "white") + 
  stat_compare_means(comparisons = list(c("D14_to_D10", "D10_to_D5"),
                                        c("D14_to_D5", "D10_to_D5"), 
                                        c("D14_to_D5", "D14_to_D10")),
                     method.args = list(alternative = "greater"), paired = TRUE) +
  theme_bw(base_size = 20) +
  rremove("legend")

ggarrange(violin1, violin2, labels = c("A", "B"))

# Test equality of IRS variance at each locus between time points ---------
Track_set_df_test <- IES_tab1 %>%
  select("IES_ID", ends_with("D5"), ends_with("D10"), ends_with("D14"))

## Bonett-Seier Test for Equality of Variability Measures
variance_test_df <- do.call(rbind, apply(Track_set_df_test, 1, 
                                         function(x) as.data.frame(Bonett.Seier.test(x = as.numeric(x[2:5]), 
                                                                                     y = as.numeric(x[10:12]), # 6:9
                                                                                     alternative = "less"))))
## bind dfs
Track_set_df_test <- as_tibble(cbind(Track_set_df_test, variance_test_df))

## inspect
Track_set_df_test %>% 
  select(-ends_with("D10"), -IES_ID) %>% 
  filter(Estimate != "NaN", Estimate != "Inf", p.value < 0.05) %>% 
  print(width = Inf, n = Inf)

# Simulation Haploid model ------------------------------------------------
## Simulate scores with haploid model
## D5 -> D14 (9 days)
## D5 -> D10 (5 days)

### compute mean at D5 (IRS0)
Track_set_df_red <- IES_tab1 %>%
  select("IES_ID", ends_with("D5"), ends_with("D14")) %>%
  mutate(IRS_P25_sc_D5 = rowMeans(across(ends_with("D5")))) %>%
  select(IES_ID, IRS_P25_sc_D5, everything()) %>%
  select(-c(3:6))

### set up dataframes
hetero_df <- as.data.frame(matrix(nrow = 0, ncol = 6))
colnames(hetero_df) <- c("GEN","IRS0", "Homoz", "H", "sd", "P_dist")
hetero_df_all <- hetero_df

### Define parameters
k = 860 # number of elements drawn (860)
Chr = 1 # number of somatic chromosomes
N = 2*k*Chr # total number of elements to draw from, number of elements after phase S (2x860)
gr = 3.476 # div. / 24h
days = 9 # re-isolation days
GEN <- seq(1, round((days * gr)), 1) # set number of generations, D9 (GEN=200, full clonal cycle)
IRS0 <- Track_set_df_red$IRS_P25_sc_D5 # Vector of IRS0 to calculate the rate of SA and H loss
#IRS0 <- 0.5
spread = "" # standard deviation

### execute senes core module
hetero_df_all <- senes_core(k, Chr, N, IRS0, GEN, hetero_df, hetero_df_all) %>%
  as_tibble() %>% 
  filter(GEN == 31)

Track_set_df_red <- bind_cols(Track_set_df_red, select(hetero_df_all, sd, P_dist)) %>%
  rename(sd_haplo = sd, P_dist_haplo = P_dist)

# Build 95% Confidence Interval haploid -----------------------------------
## Apply to all, use mapply to pass two parameters (x and SA) to LB function
Track_set_df_red$CI_95_LB_haplo <- mapply(LB, x=Track_set_df_red$IRS_P25_sc_D5, 
                                          SD=Track_set_df_red$sd_haplo)
Track_set_df_red$CI_95_UB_haplo <- mapply(UB, x=Track_set_df_red$IRS_P25_sc_D5, 
                                          SD=Track_set_df_red$sd_haplo)

### Evaluate if all IRS REPs fall within the simulated C.I.
### round up to 4 digits to perform logical evaluation
test_df <- bind_cols(select(Track_set_df_red, "IES_ID"),
                     round(select(Track_set_df_red, ends_with("D14"), matches("B")), digits = 4))

### define function at_least_one outside C.I.
at_least_one_out <- function(tab) {
  if(is.na(tab[1])){return (NA)}
  if ((tab[2] >= tab[5]) & (tab[2] <= tab[6]) &
      (tab[3] >= tab[5]) & (tab[3] <= tab[6]) &
      (tab[4] >= tab[5]) & (tab[4] <= tab[6])) {
    return(TRUE);
  }
  else{
    return(FALSE);
  }
} 

## apply at_least_one
test_df <- test_df %>%
  mutate(Test = apply(test_df, 1, at_least_one_out))

Track_set_df_red$Test <- test_df$Test

# fraction of the 75 loci for which at_least_one_out is TRUE
length(which(test_df$Test)) / nrow(test_df)

# Evaluation if at least one IRS REP falls within the simulated C.I.
# define function inside C.I.
inside_ci <- function(tab) {
  if(is.na(tab[1])){return (NA)}
  if ((tab[6] >= tab[2]) & (tab[6] <= tab[3])) {
    return(TRUE);
  }
  else{
    return(FALSE);
  }
}

# Reshape to long format. Gather REPs under a single variable
test_df <- pivot_longer(data = test_df, 
                        cols = ends_with("D14"), 
                        values_to = "IRS_P25",
                        names_to = "REP",
                        names_prefix = "IRS_P25_")

# apply inside_ci
test_df <- test_df %>%
  mutate(Test1 = apply(test_df, 1, inside_ci))

# fraction of the 75 x 3 values for which inside_ci is TRUE
length(which(test_df$Test1)) / nrow(test_df)

# Plot-1 Confidence Interval plot haploid ---------------------------------
## Modified Cleveland's dot plot
g <- ggdotchart(Track_set_df_red, x = "IES_ID", y = "IRS_P25_sc_D5",
                combine = F,
                color = "Test",                               # Color by groups
                palette = c("#FC4E07", "#00AFBB"),            # Custom color palette
                sorting = "ascending",                        # Sort value in descending order
                rotate = TRUE,                                # Rotate vertically
                dot.size = 2,                                 # Large dot size
                y.text.col = TRUE,                            # Color y text by groups
                ggtheme = theme_pubr(),                       # ggplot2 theme
                ylab = "retention levels (IRS)",
                font.label = list(size = 100)
) +
  theme_cleveland() +                                         # Add dashed grids
  theme(axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 20, vjust=-0.5),
        axis.title.x = element_text(size = 30, vjust=-2),
        plot.margin = unit(c(1,1,1,1), "cm"))                                      


g$layers <- c(geom_pointrange(mapping = aes(x = reorder(IES_ID, - IRS_P25_sc_D5), 
                                            y = IRS_P25_sc_D5, ymin = CI_95_LB_haplo, 
                                            ymax = CI_95_UB_haplo, 
                                            fill = "CI95"), 
                              data = Track_set_df_red),
              geom_point(mapping = aes(x = IES_ID, 
                                       y = IRS_P25_sc1_D14, 
                                       fill = "IRS_P25_sc1_D14"), 
                         data = Track_set_df_red,
                         fill = "orange",
                         shape = 22),
              geom_point(mapping = aes(x = IES_ID, 
                                       y = IRS_P25_sc3_D14, 
                                       fill = "IRS_P25_sc3_D14"), 
                         data = Track_set_df_red,
                         fill = "orange",
                         shape = 22),
              geom_point(mapping = aes(x = IES_ID, 
                                       y = IRS_P25_sc4_D14, 
                                       fill = "IRS_P25_sc4_D14"), 
                         data = Track_set_df_red,
                         fill = "orange",
                         shape = 22),
              g$layers)

CI_haplo <- ggpar(g, legend = "top", font.legend = c(20, "plain", "black"), 
                  legend.title = list(color = "Inside", "CI95"),
                  yticks.by = 0.1) + theme_cleveland()

# Simulation Chromosomal model --------------------------------------------
## Simulate scores with chromosomal model
# - D5 -> D14

## set up dataframes
hetero_df <- as.data.frame(matrix(nrow = 0, ncol = 6))
colnames(hetero_df) <- c("GEN","IRS0", "Homoz", "H", "sd", "P_dist")
hetero_df_all <- hetero_df

## Define parameters
k = 860 # number of elements drawn (860)
Chr = 115 # number of somatic chromosomes
N = 2*k*Chr # total number of elements to draw from, number of elements after phase S (2x860)
gr = 3.476 # div. / 24h
days = 9 # re-isolation days
GEN <- seq(1, round((days * gr)), 1) # set number of generations, D9 (GEN = ~31)
IRS0 <- Track_set_df_red$IRS_P25_sc_D5 # Vector of IRS0 to calculate the rate of SA and H loss
#IRS0 <- 0.5
spread = "" # standard deviation

## execute SENES core module
#hetero_df_all <- senes_core(k, Chr, N, IRS0, GEN, hetero_df, hetero_df_all) %>%
  #filter(GEN == 31)

#Track_set_df_red <- bind_cols(Track_set_df_red, select(hetero_df_all, sd, P_dist)) %>%
  #rename(sd_chromo = sd, P_dist_chromo = P_dist)

# Load pre-calculated estimates from file ---------------------------------
## sd_chromo
Track_set_df_red_full <- read_csv("./infiles/track_set_df_red_chromo_haplo_no_dist")
Track_set_df_red$sd_chromo <-  Track_set_df_red$sd_haplo * Track_set_df_red_full$sd_model_ratio

## Build 95% Confidence Interval chromosomal -------------------------------
### Apply to all, use mapply to pass two parameters (x and SA) to LB function
Track_set_df_red$CI_95_LB_chromo <- mapply(LB, x=Track_set_df_red$IRS_P25_sc_D5, 
                                          SD=Track_set_df_red$sd_chromo)
Track_set_df_red$CI_95_UB_chromo <- mapply(UB, x=Track_set_df_red$IRS_P25_sc_D5, 
                                          SD=Track_set_df_red$sd_chromo)

## Evaluate if all IRS REPs fall within the simulated C.I.
### round up to 4 digits to perform logical evaluation
test_df <- bind_cols(select(Track_set_df_red, "IES_ID"),
                     round(select(Track_set_df_red, ends_with("D14"), matches("B"), -matches("haplo")), digits = 4))

## apply at_least_one
test_df <- test_df %>%
  mutate(Test = apply(test_df, 1, at_least_one_out))

Track_set_df_red$Test_chromo <- test_df$Test

## fraction of the 75 loci for which at_least_one_out is TRUE
length(which(test_df$Test)) / nrow(test_df)

## Evaluation if at least one IRS REP falls within the simulated C.I.

## Reshape to long format. Gather REPs under a single variable
test_df <- pivot_longer(data = test_df, 
                        cols = ends_with("D14"), 
                        values_to = "IRS_P25",
                        names_to = "REP",
                        names_prefix = "IRS_P25_")

## apply inside_ci
test_df <- test_df %>%
  mutate(Test1 = apply(test_df, 1, inside_ci))

## fraction of the 75 x 3 values for which inside_ci is TRUE
length(which(test_df$Test1)) / nrow(test_df)

# Plot-2 Confidence Interval plot chromosomal -----------------------------
## Modified Cleveland's dot plot
g <- ggdotchart(Track_set_df_red, x = "IES_ID", y = "IRS_P25_sc_D5",
                combine = F,
                color = "Test_chromo",                               # Color by groups
                palette = c("#FC4E07", "#00AFBB"),            # Custom color palette
                sorting = "ascending",                        # Sort value in descending order
                rotate = TRUE,                                # Rotate vertically
                dot.size = 2,                                 # Large dot size
                y.text.col = TRUE,                            # Color y text by groups
                ggtheme = theme_pubr(),                       # ggplot2 theme
                ylab = "retention levels (IRS)",
                font.label = list(size = 100)
) +
  theme_cleveland() +                                         # Add dashed grids
  theme(axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 20, vjust=-0.5),
        axis.title.x = element_text(size = 30, vjust=-2),
        plot.margin = unit(c(1,1,1,1), "cm"))                                      


g$layers <- c(geom_pointrange(mapping = aes(x = reorder(IES_ID, - IRS_P25_sc_D5), 
                                            y = IRS_P25_sc_D5, ymin = CI_95_LB_chromo, 
                                            ymax = CI_95_UB_chromo, 
                                            fill = "CI95"), 
                              data = Track_set_df_red),
              geom_point(mapping = aes(x = IES_ID, 
                                       y = IRS_P25_sc1_D14, 
                                       fill = "IRS_P25_sc1_D14"), 
                         data = Track_set_df_red,
                         fill = "orange",
                         shape = 22),
              geom_point(mapping = aes(x = IES_ID, 
                                       y = IRS_P25_sc3_D14, 
                                       fill = "IRS_P25_sc3_D14"), 
                         data = Track_set_df_red,
                         fill = "orange",
                         shape = 22),
              geom_point(mapping = aes(x = IES_ID, 
                                       y = IRS_P25_sc4_D14, 
                                       fill = "IRS_P25_sc4_D14"), 
                         data = Track_set_df_red,
                         fill = "orange",
                         shape = 22),
              g$layers)

CI_chromo <- ggpar(g, legend = "top", font.legend = c(20, "plain", "black"), 
                   legend.title = list(color = "Inside", "CI95"),
                   yticks.by = 0.1) + theme_cleveland()
