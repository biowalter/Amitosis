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
library("tidyverse")
library("patchwork")

# Load data ---------------------------------------------------------------
COV_W_combo_All <- read.table("./infiles/COV_W_combo_All.txt", 
                                header = T, stringsAsFactors = F, sep = "\t")

# Remove internal headers after concatenation
indi <- as.numeric(rownames(COV_W_combo_All[COV_W_combo_All$W_ID=="W_ID", ]))
COV_W_combo_All <- COV_W_combo_All[-indi, ]
# Remove failed sc sample
COV_W_combo_All <- COV_W_combo_All[-grep("sc2_D14", COV_W_combo_All$sample), ]

## NGS stats
ngs <- read_tsv("./infiles/NGS_stats_header")

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

S2B <- ggplot() + 
  scale_y_continuous(name = "coverage", breaks = seq(0, 120, 10), limits = c(0, 120)) + # add secondary axis
  geom_boxplot(mapping = aes(x=W_ID, y=W_COV, fill = "#008080"), data = COV_W_combo_All[COV_W_combo_All$sample!="aDNA" & COV_W_combo_All$sample!="mcDNA", ], outlier.color = "gray") + # box plot for scDNA samples
  geom_boxplot(mapping = aes(x=W_ID, y=W_COV, fill = "dark orange"), data = COV_W_combo_All[COV_W_combo_All$sample=="mcDNA", ], alpha = 0.8, colour="black", outlier.fill = "dark orange") +
  geom_boxplot(mapping = aes(x=W_ID, y=W_COV, fill = "sky blue"), data = COV_W_combo_All[COV_W_combo_All$sample=="aDNA", ], alpha = 1.0, colour = "navy blue") +
  scale_fill_identity(name = "Samples",
                      breaks = c("#008080", "sky blue", "dark orange"),
                      labels = c("scDNA", "aDNA", "mcDNA"),
                      guide = "legend") +
  theme_bw(base_size = 22) +
  scale_x_discrete(name = "distance from chromosome ends (kb)", breaks = as.character(seq(0, 50, 5)), limits = as.character(seq(1, 30, 1))) +
  theme(legend.position="top")

S2B

# Effect of increase COV on under-represented regions ---------------------
## add average COV information to data frame
ngs <- ngs %>% 
  select(Sample, COV, `Mapped pairs (M)`)

COV_dist_full <- COV_W_combo_All
sample_index <- lapply(ngs$Sample, function(x) grep(x, COV_dist_full$sample))

for (i in seq_along(ngs$Sample)) {
  COV_dist_full$AVG_COV[sample_index[[i]]] <- ngs$COV[i]
}

## plot COV against AVG_COV by DIST
### filter GC content bins
DIST_bins <- seq(5, 30, 5)
COV_dist_full_bins <- COV_dist_full %>% 
  filter(W_ID %in% DIST_bins)

### reorder GC bin levels
COV_dist_full_bins$W_ID <- factor(COV_dist_full_bins$W_ID, levels = as.character(DIST_bins))

### add reference y-intercept to annotate facets with mcDNA coverage
DIST_index <- lapply(DIST_bins, function(x) grep(paste("^", x, sep = ""), COV_dist_full_bins$W_ID))
REF_W_COV_vector <- COV_dist_full_bins %>% 
  filter(sample == "mcDNA") %>% 
  select(W_COV) %>% 
  unlist() %>% 
  as.vector()

COV_dist_full_bins$REF_COV <- NA
for (i in seq_along(DIST_index)) {
  COV_dist_full_bins$REF_COV[DIST_index[[i]]] <- REF_W_COV_vector[i]
}

## compute average by DIST-bin (W_ID)
COV_dist_full_bins %>% 
  filter(sample == "mcDNA") %>% 
  group_by(W_ID) %>% 
  summarise(REF_COV_avg = mean(REF_COV, na.rm = TRUE)) %>% 
  print(n = Inf)

### ref x-intercept 1.5 mcDNA COV
x_intercept <- 1.5 * COV_dist_full_bins %>% 
  filter(sample == "mcDNA") %>% 
  select(AVG_COV) %>% 
  unique() %>% 
  as.numeric()

### exclude reference samples
COV_dist_full_bins <- COV_dist_full_bins %>%
  filter(sample != "aDNA", sample != "mcDNA")

### color palette
pal <- colorRampPalette(pal_jama()(7))(length(levels(as.factor(COV_dist_full_bins$sample))))

### plot
S2D <- COV_dist_full_bins %>%
  ggplot(mapping = aes(x = AVG_COV, y = W_COV, fill = sample)) +
  geom_boxplot(width = 5, varwidth = T) +
  geom_hline(mapping = aes(yintercept = REF_COV),
             linetype = "dashed", color = "dark orange", size = 1) +
  geom_vline(xintercept = x_intercept,
             linetype = "dotted", color = "black", size = 0.8) +
  xlab("coverage (genome)") +
  ylab("coverage (kb-bin)") +
  scale_fill_manual(values = pal) +
  facet_wrap(~ W_ID, ncol = 3) +
  ylim(c(0, 100)) +
  theme_bw(base_size = 22) +
  theme(legend.position="top", 
        legend.key.size = unit(14, "points"),
        legend.text = element_text(size = 14),
        legend.spacing.x = unit(5, "points"))

S2D

### patch FigS2 in 
S2A + S2B + S2C + S2D + plot_annotation(tag_levels = 'A')

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
