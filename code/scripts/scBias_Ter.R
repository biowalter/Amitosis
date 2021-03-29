# Description -------------------------------------------------------------
# Authors: Vitali et al. (2021)
# Quantitative analysis of genome representation biases II
## Terminal bias
# .........................................................................

# Libraries ---------------------------------------------------------------
library("ggplot2")
library("ggpubr")
library("RColorBrewer")
library("ggsci")

# Load data ---------------------------------------------------------------
COV_W_combo_All <- read.table("./infiles/COV_W_combo_All.txt", 
                                header = T, stringsAsFactors = F, sep = "\t")

# Remove internal headers after concatenation
indi <- as.numeric(rownames(COV_W_combo_All[COV_W_combo_All$W_ID=="W_ID", ]))
COV_W_combo_All <- COV_W_combo_All[-indi, ]
# Remove failed sc sample
COV_W_combo_All <- COV_W_combo_All[-grep("sc2_D14", COV_W_combo_All$sample), ]

# Plots  ------------------------------------------------------------------
# - Plot-1 COV against distance from chromosome ends (W)
# - Plot-2 Norm. COV against distance from chromosome ends (W)

# Plot1  ------------------------------------------------------------------
## COV by distance from chromosome termini
# - 2kb tiles centered every kb

COV_W_combo_All$W_ID <- factor(x = COV_W_combo_All$W_ID, levels = as.character(seq(1:49)))
# Reorder levels (windows in ascending order 1 to 49)
COV_W_combo_All$W_ID <- factor(COV_W_combo_All$W_ID, levels = order(COV_W_combo_All$W_ID))
COV_W_combo_All$W_COV <-as.numeric(COV_W_combo_All$W_COV)

p <- ggplot() + 
  scale_y_continuous(name = "Coverage (reads / base)", breaks = seq(0, 120, 10), limits = c(0, 120)) + # add secondary axis
  geom_boxplot(mapping = aes(x=W_ID, y=W_COV, fill = "W_ID"), data = COV_W_combo_All[COV_W_combo_All$sample!="aDNA" & COV_W_combo_All$sample!="mcDNA", ], outlier.color = "gray") + # box plot for scDNA samples
  scale_fill_brewer(palette = "Dark2") +
  #geom_jitter(mapping = aes(x=W_ID, y=W_COV, fill = "W_ID"), data = COV_W_combo_All[COV_W_combo_All$sample!="aDNA" & COV_W_combo_All$sample!="mcDNA", ], shape=16, position=position_jitter(0.2)) +
  geom_boxplot(mapping = aes(x=W_ID, y=W_COV, fill = "W_ID"), data = COV_W_combo_All[COV_W_combo_All$sample=="aDNA", ], alpha = 1.0, colour = "navy blue", fill = "sky blue") +
  geom_boxplot(mapping = aes(x=W_ID, y=W_COV, fill = "W_ID"), data = COV_W_combo_All[COV_W_combo_All$sample=="mcDNA", ], alpha = 0.2, colour="black", fill = "dark orange", outlier.fill = "dark orange") +
  theme_bw(base_size = 22) +
  scale_x_discrete(name = "d from chr ends (kbp)", breaks = as.character(seq(0, 50, 5)), limits = as.character(seq(1, 30, 1)))
p

# Plot-2 ------------------------------------------------------------------
## Norm. COV by distance from chromosome termini
#- 2kb tiles centered every kb

COV_W_combo_All$W_ID <- factor(x = COV_W_combo_All$W_ID, levels = as.character(seq(1:49)))
# Reorder levels (windows in ascending order 1 to 49)
COV_W_combo_All$W_ID <- factor(COV_W_combo_All$W_ID, levels = order(COV_W_combo_All$W_ID))
COV_W_combo_All$COV_norm <-as.numeric(COV_W_combo_All$COV_norm)

p1 <- ggplot() + 
  scale_y_continuous(name = "normalized coverage", breaks = seq(0, 2.0, 0.2), limits = c(0, 2.0)) + # add secondary axis
  geom_boxplot(mapping = aes(x=W_ID, y=COV_norm, fill = "#00468BFF"), data = COV_W_combo_All[COV_W_combo_All$sample!="aDNA" & COV_W_combo_All$sample!="mcDNA", ], outlier.color = "gray") + # box plot for scDNA samples
  #scale_fill_lancet() +
  geom_boxplot(mapping = aes(x=W_ID, y=COV_norm, fill = "sky blue"), data = COV_W_combo_All[COV_W_combo_All$sample=="aDNA", ], alpha = 1.0, colour = "navy blue") +
  geom_boxplot(mapping = aes(x=W_ID, y=COV_norm, fill = "dark orange"), data = COV_W_combo_All[COV_W_combo_All$sample=="mcDNA", ], alpha = 0.2, colour="black", outlier.fill = "dark orange") +
  scale_fill_identity(name = "Samples",
                       breaks = c("#00468BFF", "sky blue", "dark orange"),
                       labels = c("scDNA", "aDNA", "mcDNA"),
                       guide = "legend") +
  theme_bw(base_size = 25) +
  theme(legend.position="top", legend.box = "horizontal",
        plot.margin = unit(c(1,3,1,1), "cm")) +
  scale_x_discrete(name = "distance from chromosome ends (kb)", breaks = as.character(seq(0, 50, 5)), limits = as.character(seq(1, 30, 1)))
Fig1b <- p1

# Terminal bias -----------------------------------------------------------
# Genome representation bias 2. 
# - Term. bias-1 = COV_30kb / COV_1kb
# - Term. bias-2 = inear regression of norm. coverage vs. distance (slope)

# Term. bias-1 ------------------------------------------------------------
## COV_30kb / COV_1kb

COV_W_combo_All$W_ID <- as.numeric(as.character(COV_W_combo_All$W_ID))
## 30 kb away from chr termini
above <- mean(COV_W_combo_All$COV_norm[COV_W_combo_All$W_ID==30])
## 1kb kb away from chr termini
below <- mean(COV_W_combo_All$COV_norm[COV_W_combo_All$W_ID==1])
## W_ID-bias metric
Ter_bias <- above / below
Ter_bias
## Compute for each sample individually
SAMP <- levels(as.factor(COV_W_combo_All$sample))
Ter_bias_df <- as.data.frame(matrix(NA, 1, 4))
colnames(Ter_bias_df) <- c("Ter_bias", "COV_1kb_obs", "COV_30kb_obs", "sample")
for (i in 1:length(SAMP)) {
  COV_W_combo_All1 <- COV_W_combo_All[COV_W_combo_All$sample==SAMP[i], ]
  rownames(COV_W_combo_All1[COV_W_combo_All1$sample==SAMP[i], ]) <- seq(1, 11270, 1)
  above <- mean(COV_W_combo_All1$COV_norm[COV_W_combo_All1$W_ID==30])
  below <- mean(COV_W_combo_All1$COV_norm[COV_W_combo_All1$W_ID==1])
  Ter_bias <- above / below
  Ter_bias_df[i, ] <- cbind(Ter_bias, below, above, SAMP[i])
}

Ter_bias_df
Ter_bias_df[, 1:3] <- apply(Ter_bias_df[, 1:3], 2, as.numeric)

# Term. bias-2 ------------------------------------------------------------
## slope of regression (Norm_COV on distance)

COV_W_combo_All$W_ID <- COV_W_combo_All$W_ID / 10 # Norm_COV every 10kb
summary(COV_W_combo_All$W_ID)
# Linear regression for aDNA
aDNA_corr <- COV_W_combo_All[grep("aDNA", COV_W_combo_All$sample),]
lm(aDNA_corr$COV_norm[aDNA_corr$W_ID > 0.1 & aDNA_corr$W_ID < 3] ~ aDNA_corr$W_ID[aDNA_corr$W_ID > 0.1 & aDNA_corr$W_ID < 3])
# Correlation in mcDNA
mcDNA_corr <- COV_W_combo_All[grep("mcDNA", COV_W_combo_All$sample),]
lm(mcDNA_corr$COV_norm[mcDNA_corr$W_ID > 0.1 & mcDNA_corr$W_ID < 3] ~ mcDNA_corr$W_ID[mcDNA_corr$W_ID > 0.1 & mcDNA_corr$W_ID < 3])
# Compute for all
lm(COV_W_combo_All$COV_norm[COV_W_combo_All > 0.1 & COV_W_combo_All$W_ID < 3] ~ COV_W_combo_All$W_ID[COV_W_combo_All > 0.1 & COV_W_combo_All$W_ID < 3])
## Compute for each sample individually
SAMP <- levels(as.factor(COV_W_combo_All$sample))
Ter_bias_df1 <- as.data.frame(matrix(NA, 1, 4))
colnames(Ter_bias_df1) <- c("Ter_bias", "COV_1kb_est", "COV_30kb_est", "sample")
for (i in 1:length(SAMP)) {
  COV_W_combo_All1 <- COV_W_combo_All[COV_W_combo_All$sample==SAMP[i], ]
  rownames(COV_W_combo_All1[COV_W_combo_All1$sample==SAMP[i], ]) <- seq(1, 11270, 1)
  Ter_bias <- lm(COV_W_combo_All1$COV_norm[COV_W_combo_All1 > 0.1 & COV_W_combo_All1$W_ID < 3] ~ COV_W_combo_All1$W_ID[COV_W_combo_All1 > 0.1 & COV_W_combo_All1$W_ID < 3])[[1]][2] # Slope of lm COV on W_ID%
  intercept <- lm(COV_W_combo_All1$COV_norm[COV_W_combo_All1 > 0.1 & COV_W_combo_All1$W_ID < 3] ~ COV_W_combo_All1$W_ID[COV_W_combo_All1 > 0.1 & COV_W_combo_All1$W_ID < 3])[[1]][1]
  above <- (Ter_bias * 3) + intercept
  below <- (Ter_bias * 0.1) + intercept
  Ter_bias_df1[i, ] <- cbind(Ter_bias, below, above, SAMP[i])
}

Ter_bias_df1
Ter_bias_df1[, 1:3] <- apply(Ter_bias_df1[, 1:3], 2, as.numeric)

# Estimated vs observed bias ----------------------------------------------
## compare observed bias with estimates from linear regression

Ter_bias_combo <- cbind(Ter_bias_df1, Ter_bias_df[, 2:3])
