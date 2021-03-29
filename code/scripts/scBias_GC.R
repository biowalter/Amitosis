# Description -------------------------------------------------------------
# Authors: Vitali et al. (2021)
# Quantitative analysis of genome representation biases I
## GC bias
# .........................................................................

# Libraries ---------------------------------------------------------------
library("ggplot2")
library("ggpubr")
library("RColorBrewer")
library("ggsci")
library("radiant.data")
library("dplyr")

# Load data ---------------------------------------------------------------
combo_ID_COV_SIZE <- read.table("./infiles/combo_ID_COV_SIZE.txt", 
                                header = T, stringsAsFactors = F, sep = "\t")

sc_combo <- read.table("./infiles/sc_combo.txt", header = T,
                       stringsAsFactors = F, sep = "\t", quote = "", dec = ".")

aDNA <- sc_combo[grep(paste("^", "aDNA", sep = ""), sc_combo$sample), ]
mcDNA <- sc_combo[grep(paste("^", "mcDNA", sep = ""), sc_combo$sample), ]
ctr_df <- rbind(aDNA, mcDNA)
sc_combo <- sc_combo[grep(paste("^", "sc", sep = ""), sc_combo$sample), ]


# Plots -------------------------------------------------------------------
# - Plot-1 coverage (COV) against GC content (GC)
# - Plot-2 Normalized coverage (norm COV) against GC content (GC)

# Plot-1 ------------------------------------------------------------------
## Coverage (COV) against GC content (GC)

sc_combo <- sc_combo[-grep("sc2_D14", sc_combo$sample),]
sc_combo$GC <- as.factor(sc_combo$GC)
aDNA$GC <- as.factor(aDNA$GC)
mcDNA$GC <- as.factor(mcDNA$GC) 
p <- ggplot() + 
  scale_y_continuous(name = "Coverage (reads / base)", breaks = seq(0, 120, 10), sec.axis = sec_axis(~ . * (40000) / (10^6)  , name = " # Mb DNA (GC-bin)"), limits = c(0, 120)) + # add secondary axis
  geom_boxplot(mapping = aes(x=GC, y=COV, fill = "GC"), data = sc_combo) + # box plot for scDNA samples
  scale_fill_brewer(palette = "Dark2") +
  geom_jitter(mapping = aes(x=GC, y=COV, fill = "GC"), data = sc_combo, shape=16, position=position_jitter(0.2)) +
  geom_col(mapping = aes(x=GC, y=WINDOWS/40000), data = aDNA, alpha = 0.1, color="black", fill="white") + # GC distribution bar chart 
  geom_line(mapping = aes(x=GC, y=COV, group = "GC"), data = aDNA, alpha = 0.8, color="dark blue", size = 1) +
  geom_line(mapping = aes(x=GC, y=COV, group = "GC"), data = mcDNA, alpha = 0.8, color="dark orange", size = 1) +
  theme_bw(base_size = 22) +
  scale_x_discrete(name = "GC content (%)", breaks = as.character(seq(0, 55, 5)), limits = as.character(seq(0, 55, 1)))
p

# Plot-2 ------------------------------------------------------------------
## Normalized coverage (norm COV) against GC content (GC)

sc_combo$GC <- as.factor(sc_combo$GC)
aDNA$GC <- as.factor(aDNA$GC)
mcDNA$GC <- as.factor(mcDNA$GC) 
p <- ggplot() + 
  scale_y_continuous(name = "Normalized COV (reads / base / Mean COV)", breaks = seq(0, 2.0, 0.2), sec.axis = sec_axis(~ . * 2500000 / (10^6)  , name = " # Mb DNA (GC-bin)"), limits = c(0, 2.0)) + # add secondary axis
  geom_boxplot(mapping = aes(x=GC, y=NORMALIZED_COVERAGE, fill = "#00468BFF"), data = sc_combo) + # box plot for scDNA samples
  scale_fill_identity(name = "scDNA",
                       breaks = c("#00468BFF"),
                       labels = c("all (n=11)"),
                       guide = "legend") +
  geom_jitter(mapping = aes(x=GC, y=NORMALIZED_COVERAGE, fill = "#00468BFF"), data = sc_combo, shape=16, position=position_jitter(0.2)) +
  geom_col(mapping = aes(x=GC, y=WINDOWS/2500000), data = aDNA, alpha = 0.1, color="black", fill="white") + # GC distribution bar chart 
  geom_line(mapping = aes(x=GC, y=NORMALIZED_COVERAGE, group = "GC", color="dark blue"), data = aDNA, alpha = 0.8, size = 1) +
  geom_line(mapping = aes(x=GC, y=NORMALIZED_COVERAGE, group = "GC", color="dark orange"), data = mcDNA, alpha = 0.8, size = 1) +
  scale_color_identity(name = "Control",
                       breaks = c("dark blue", "dark orange"),
                       labels = c("aDNA", "mcDNA"),
                       guide = "legend") +
  theme_bw(base_size = 22) +
  scale_x_discrete(name = "GC content (%)", breaks = as.character(seq(0, 55, 5)), limits = as.character(seq(0, 55, 1)))

p

# GC bias -----------------------------------------------------------------
# - GC bias-1 = linear regression of coverage vs. GC (slope)
# - GC bias-2 = (COV_mean + sd) / (COV_mean - sd)
# - GC bias-3 = linear regression of norm. coverage vs. GC (slope)

# GC bias-1 ---------------------------------------------------------------
## slope of regression (coverage)

sc_combo <- rbind(ctr_df, sc_combo)
sc_combo$GC <- as.numeric(as.character(sc_combo$GC))
# mean percent GC content 
w_mean <- weighted.mean(x = sc_combo$GC[sc_combo$sample=="mcDNA"], w = sc_combo$WINDOWS[sc_combo$sample=="mcDNA"])
# sd percent GC content 
w_sd <- weighted.sd(sc_combo$GC[sc_combo$sample=="mcDNA"], w = sc_combo$WINDOWS[sc_combo$sample=="mcDNA"])
# Linear regression for aDNA
aDNA_corr <- sc_combo[grep("aDNA", sc_combo$sample),]
lm(aDNA_corr$COV[aDNA_corr$GC > 9 & aDNA_corr$GC < 50] ~ aDNA_corr$GC[aDNA_corr$GC > 9 & aDNA_corr$GC < 50])
# Correlation in mcDNA
mcDNA_corr <- sc_combo[grep("mcDNA", sc_combo$sample),]
lm(mcDNA_corr$COV[mcDNA_corr$GC > 9 & mcDNA_corr$GC < 50] ~ mcDNA_corr$GC[mcDNA_corr$GC > 9 & mcDNA_corr$GC < 50])
# Compute for all
lm(sc_combo$COV[sc_combo$GC > 9 & sc_combo$GC < 50] ~ sc_combo$GC[sc_combo$GC > 9 & sc_combo$GC < 50])
## Compute for each sample individually
SAMP <- levels(as.factor(sc_combo$sample))
GC_bias_df <- as.data.frame(matrix(NA, 1, 4))
colnames(GC_bias_df) <- c("GC_bias", "COV_mu-sd.", "COV_mu+sd", "sample")
for (i in 1:length(SAMP)) {
  sc_combo1 <- sc_combo[sc_combo$sample==SAMP[i], ]
  rownames(sc_combo1[sc_combo1$sample==SAMP[i], ]) <- seq(1, 101, 1)
  GC_bias <- lm(sc_combo1$COV[sc_combo1$GC > 9 & sc_combo1$GC < 50] ~ sc_combo1$GC[sc_combo1$GC > 9 & sc_combo1$GC < 50])[[1]][2] # Slope of lm COV on GC%
  print(GC_bias)
  intercept <- lm(sc_combo1$COV[sc_combo1$GC > 9 & sc_combo1$GC < 50] ~ sc_combo1$GC[sc_combo1$GC > 9 & sc_combo1$GC < 50])[[1]][1]
  print(intercept)
  above <- (GC_bias * (w_mean + 1 * w_sd)) + intercept
  below <- (GC_bias * (w_mean - 1 * w_sd)) + intercept
  GC_bias_df[i, ] <- cbind(GC_bias, below, above, SAMP[i])
}

GC_bias_df

# GC bias-2 ---------------------------------------------------------------
##  (COV_mean + sd) / (COV_mean - sd)

sc_combo$GC <- sc_combo$GC / 10 # Norm_COV every 10% GC
sc_combo$GC <- as.numeric(as.character(sc_combo$GC))
## 1st Qu.
above <- mean(sc_combo$NORMALIZED_COVERAGE[near(x = sc_combo$GC, y =  (w_mean + 1 * w_sd) / 10, tol = 0.1)])
## 3rd Qu.
below <- mean(sc_combo$NORMALIZED_COVERAGE[near(x = sc_combo$GC, y =  (w_mean - 1 * w_sd) / 10, tol = 0.1)])
## GC-bias metric
GC_bias <- above / below
GC_bias
## Compute for each sample individually
SAMP <- levels(as.factor(sc_combo$sample))
GC_bias_df <- as.data.frame(matrix(NA, 1, 4))
colnames(GC_bias_df) <- c("GC_bias", "COV_mu-sd_obs", "COV_mu+sd_obs", "sample")
for (i in 1:length(SAMP)) {
  sc_combo1 <- sc_combo[sc_combo$sample==SAMP[i], ]
  rownames(sc_combo1[sc_combo1$sample==SAMP[i], ]) <- seq(1, 101, 1)
  above[i] <- sc_combo1$NORMALIZED_COVERAGE[near(x = sc_combo$GC, y =  (w_mean + 1 * w_sd) / 10, tol = 0.1)]
  below[i] <- sc_combo1$NORMALIZED_COVERAGE[near(x = sc_combo$GC, y =  (w_mean - 1 * w_sd) / 10, tol = 0.1)]
  GC_bias[i] <- above[i] / below[i]
  GC_bias_df[i, ] <- cbind(GC_bias[i], below[i], above[i], SAMP[i])
}

GC_bias_df
GC_bias_df[, 1:3] <- apply(GC_bias_df[, 1:3], 2, as.numeric)


# GC bias-3 ---------------------------------------------------------------
##  slope of regression (normalized coverage)

# sc_combo$GC <- sc_combo$GC / 10 # Norm_COV every 10% GC
# Linear regression for aDNA
aDNA_corr <- sc_combo[grep("aDNA", sc_combo$sample),]
lm(aDNA_corr$NORMALIZED_COVERAGE[aDNA_corr$GC > 0.9 & aDNA_corr$GC < 5] ~ aDNA_corr$GC[aDNA_corr$GC > 0.9 & aDNA_corr$GC < 5])
# Correlation in mcDNA
mcDNA_corr <- sc_combo[grep("mcDNA", sc_combo$sample),]
lm(mcDNA_corr$NORMALIZED_COVERAGE[mcDNA_corr$GC > 0.9 & mcDNA_corr$GC < 5] ~ mcDNA_corr$GC[mcDNA_corr$GC > 0.9 & mcDNA_corr$GC < 5])
# Compute for all
lm(sc_combo$NORMALIZED_COVERAGE[sc_combo$GC > 0.9 & sc_combo$GC < 5] ~ sc_combo$GC[sc_combo$GC > 0.9 & sc_combo$GC < 5])
## Compute for each sample individually
SAMP <- levels(as.factor(sc_combo$sample))
GC_bias_df1 <- as.data.frame(matrix(NA, 1, 4))
colnames(GC_bias_df1) <- c("GC_bias", "COV_mu-sd_est", "COV_mu+sd_est", "sample")
for (i in 1:length(SAMP)) {
  sc_combo1 <- sc_combo[sc_combo$sample==SAMP[i], ]
  rownames(sc_combo1[sc_combo1$sample==SAMP[i], ]) <- seq(1, 101, 1)
  model = lm(sc_combo1$NORMALIZED_COVERAGE[sc_combo1$GC > 0.9 & sc_combo1$GC < 5] ~ sc_combo1$GC[sc_combo1$GC > 0.9 & sc_combo1$GC < 5])[[1]]
  GC_bias <- model[2] # Slope of COV against GC%
  intercept <- model[1] # intercept
  above <- (GC_bias * (w_mean + 1 * w_sd) / 10) + intercept
  below <- (GC_bias * (w_mean - 1 * w_sd) / 10) + intercept
  GC_bias_df1[i, ] <- cbind(GC_bias, below, above, SAMP[i])
}

GC_bias_df1
GC_bias_df1[, 1:3] <- apply(GC_bias_df1[, 1:3], 2, as.numeric)

# Estimated vs observed bias ----------------------------------------------
## compare observed bias with estimates from linear regression

GC_bias_combo <- cbind(GC_bias_df1, GC_bias_df[, 2:3])
GC_bias_combo$sample[3:nrow(GC_bias_combo)] <- rep("sc", 11)

# calculate summary stats by group
GC_bias_combo %>% 
  group_by(sample) %>% 
  summarise(across(everything(), c(mean, sd))) %>% 
  print(width = Inf)
