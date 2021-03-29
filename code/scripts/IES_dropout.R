# Description -------------------------------------------------------------
# Authors: Vitali et al. (2021)
# Quantitative Analysis of IES dropout
# .........................................................................

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(ggplot2)

# Import data --------------------------------------------------------------
## IES table
IES_tab <- read.table("./infiles/amitosis_IES_table.txt", 
                      header = T, stringsAsFactors = F, sep = "\t", quote = "")
## Scaffolds info
scaffold_df <- read.table("./infiles/combo_ID_COV_SIZE.txt",
                          header = T, stringsAsFactors = F, sep = "\t", quote = "")
## NGS stats
NGS_stats <- read.table("./infiles/NGS_stats.txt",
                        header = T, stringsAsFactors = F, sep = "\t", quote = "")
### Rm aDNA sample (artifial reads)
NGS_stats <- NGS_stats[-1, ]

# Retrieve Terminal IESs --------------------------------------------------
## Build frame (IES_ID; IES_coor; scaffold)
IES_ter_df <- data.frame(IES_ID = IES_tab$IES_ID[!is.na(IES_tab$IES_ID) & !duplicated(IES_tab$IES_ID)],
                         IES_coor = unlist(lapply(IES_tab$IES_ID[!is.na(IES_tab$IES_ID) & !duplicated(IES_tab$IES_ID)], 
                                                  function(x) unlist(strsplit(x, split = ".", fixed = T))[5])),
                         scaffold = unlist(lapply(IES_tab$IES_ID[!is.na(IES_tab$IES_ID) & !duplicated(IES_tab$IES_ID)], 
                                                  function(x) unlist(strsplit(x, split = ".", fixed = T))[4])),
                         stringsAsFactors = F)

## Convert IES_coor to numeric
IES_ter_df$IES_coor <- as.numeric(IES_ter_df$IES_coor)

## 1. Select IESs within 30kb from chr start
IES_ter_df_sx <- IES_ter_df[IES_ter_df$IES_coor <= 30*10^3, ]

## Fetch scaffold size
ID_SIZE_df <- scaffold_df[, 1:2]
ID_SIZE_df$scaffold <- unlist(lapply(ID_SIZE_df[, 1], 
                                     function(x) unlist(strsplit(x, split = "_", fixed = T))[2]))

IES_ter_df$size_bp <- ID_SIZE_df$scaffold_size[match(IES_ter_df$scaffold, ID_SIZE_df$scaffold)]

## 2. select IESs within 30kb from chr end
IES_ter_df_rx <- IES_ter_df[(as.numeric(IES_ter_df$size_bp) - as.numeric(IES_ter_df$IES_coor)) <= 30*10^3, ]

## 2. combine sx and rx
IES_ter_df <- rbind(IES_ter_df_sx, IES_ter_df_rx[1:3])

# Build dropout frame -----------------------------------------------------
## Sample; Total; Terminal
## Total IES dropout = fraction of all IES loci (n=44,928) with coverage <= 20
## Terminal IES dropout = fraction of IES loci within 30kb of either scaffold ends
## with coverage <= 20
# .........................................................................
drop_out_df <- data.frame(sample = unique(scaffold_df$sample)[-1],
                          Total_IES_drop_out = NA,
                          Ter_IES_drop_out = NA,
                          stringsAsFactors = FALSE)

## Relabel mcDNA -> 2b_D7
drop_out_df$sample[1] <- "2b_D7"

## Coverage field vector
MINUS_COV <- unlist(lapply(drop_out_df$sample, function(x) paste("IES_MINUS_P25", x, sep = "_")))

## IES loci count (all)
IES_all <- length(IES_tab$IES_ID[!is.na(IES_tab$IES_ID) & !duplicated(IES_tab$IES_ID)])

## Compute Total & Terminal dropout
for (i in seq_along(MINUS_COV)) {
  print(grep(MINUS_COV[i], names(IES_tab), value = T))
  drop_out_df$Total_IES_drop_out[i] <- 1 - (length(IES_tab$IES_ID[!is.na(IES_tab$IES_ID) & !duplicated(IES_tab$IES_ID) & 
                                                                    IES_tab[grep(MINUS_COV[i], names(IES_tab))] > 20]) / IES_all)
  
  drop_out_df$Ter_IES_drop_out[i] <- length(which(ifelse(IES_tab[grep(MINUS_COV[i], names(IES_tab))][, 1][which(IES_tab$IES_ID%in%IES_ter_df$IES_ID)] <= 20, TRUE, FALSE))) / IES_all
}

## Rm failed scDNA sample (sc2_D14) and print frame
(drop_out_df <- drop_out_df[-grep("sc2_D14", drop_out_df$sample), ])

## Relabel 2b_D7 -> mcDNA
drop_out_df$sample[1] <- "mcDNA"

# Residual IES dropout --------------------------------------------------------
## Compute residual drop out for mcDNA, possibly unrelated to amplification biases
mcDNA_residual <- drop_out_df$Total_IES_drop_out[1] - drop_out_df$Ter_IES_drop_out[1]

## scDNA Residual dropout:
## mcDNA residual dropout scaled on the number of mapped reads (sc / mc)
mapped_reads_collate <- NGS_stats$Read_pairs_mapped[match(drop_out_df$sample, NGS_stats$sample)]
drop_out_df$Read_pairs_mapped <- mapped_reads_collate
drop_out_df$Residual_IES_drop_out <- mcDNA_residual * (mapped_reads_collate[1] / mapped_reads_collate)

# GC IES dropout ----------------------------------------------------------
## Estimate GC drop out, GC dropout = Total dropout − (Terminal + Residual)
drop_out_df$GC_IES_drop_out <- drop_out_df$Total_IES_drop_out - (drop_out_df$Ter_IES_drop_out + drop_out_df$Residual_IES_drop_out)
drop_out_df

# Quantitative analysis of IES dropout ------------------------------------
## Aggregate Stats for scDNA samples
## Compute mean
apply(drop_out_df[2:12, 2:length(drop_out_df)], 2,  mean)
## Compute sd by group
apply(drop_out_df[2:12, 2:length(drop_out_df)], 2,  sd)
## For sc with comparable number of mapped reads (between 5M and 15M)
apply(drop_out_df[drop_out_df$Read_pairs_mapped < 15*10^6 & drop_out_df$Read_pairs_mapped  > 5*10^6, ][-1, 2:6], 2, mean)
apply(drop_out_df[drop_out_df$Read_pairs_mapped < 15*10^6 & drop_out_df$Read_pairs_mapped  > 5*10^6, ][-1, 2:6], 2, sd)
## For sc with twice as many mapped reads [over 19M]
apply(drop_out_df[drop_out_df$Read_pairs_mapped > 19*10^6, ][2:6], 2, mean)
apply(drop_out_df[drop_out_df$Read_pairs_mapped > 19*10^6, ][2:6], 2, sd)

## Table 2
Tab2 <- rbind(drop_out_df[1, ], 
              c("scDNA_1x_mean", round( apply(drop_out_df[drop_out_df$Read_pairs_mapped < 15*10^6 & drop_out_df$Read_pairs_mapped  > 5*10^6, ][-1, 2:6], 2, mean), 2)),
              c("scDNA_1x_sd", round(apply(drop_out_df[drop_out_df$Read_pairs_mapped < 15*10^6 & drop_out_df$Read_pairs_mapped  > 5*10^6, ][-1, 2:6], 2, sd), 2)),
              c("scDNA_2x_mean", round(apply(drop_out_df[drop_out_df$Read_pairs_mapped > 19*10^6, ][2:6], 2, mean), 2)),
              c("scDNA_2x_sd", round(apply(drop_out_df[drop_out_df$Read_pairs_mapped > 19*10^6, ][2:6], 2, sd), 2)))

Tab2$Read_pairs_mapped <- round(as.numeric(Tab2$Read_pairs_mapped) / 10^6, 2)
Tab2

# Somatic IESs count (detected) -------------------------------------------
# somatic IESs, IRS > 0.1, IES_COV > 20

## Subset IES table
IES_tab1 <- IES_tab %>%
  filter(!is.na(IES_ID) & !duplicated(IES_ID)) %>% # rm NAs, drop duplicates
  select(starts_with("IES_ID"), matches("P25"), -contains(c("PLUS", "MINUS"))) # select relevant columns
  
## Reshape frame to long format (double pivoting on two variables)
### Pivot on IRS
pivot1 <- pivot_longer(data = IES_tab1,
                       cols = starts_with("IRS_P25"),
                       values_to = "IRS",
                       names_to = "sample",
                       names_prefix = "IRS_") %>% 
  print(n = 20, width = Inf)

### Pivot on COV
pivot2 <- pivot_longer(data = pivot1,
                       cols = starts_with("IES_COV_P25"),
                       values_to = "COV",
                       names_to = "sample1",
                       names_prefix = "IES_COV_") %>% 
  print(n = 20, width = Inf)

### Compute counts by sample
count <- pivot2 %>%
  filter(IRS > 0.1, COV > 20) %>%
  group_by(sample, sample1) %>% # group by both sample vars
  summarize(count = n())

### filter at sample == sample1
#### Add logical
count$logic <- as.vector(ifelse(count[, 1] == count[, 2], TRUE, FALSE))
#### Filter at logical
somatic <- count %>%
  filter(logic == TRUE)

## add somatic IES counts
### fix sample names
somatic$sample <- unlist(lapply(somatic$sample, function(x) paste(unlist(strsplit(x, split = "_"))[2:3], collapse = "_")))
somatic$sample[1] <- "mcDNA"

### extend dropout frame
drop_out_df_extended <- bind_cols(drop_out_df, 
                                  somatic[match(drop_out_df$sample, somatic$sample), 3])

### Change name count -> sIESs
colnames(drop_out_df_extended)[grep("count", colnames(drop_out_df_extended))] <- "sIESs"

# Fig2A -------------------------------------------------------------------
# Corr btw IES dropout and num of somatic IESs
cor.test(drop_out_df_extended$Total_IES_drop_out, as.numeric(drop_out_df_extended$sIESs))

plot(drop_out_df_extended$Total_IES_drop_out, 
     as.numeric(drop_out_df_extended$sIESs), xlim = c(0, 0.5), ylim = c(0, 400))

drop_out_df_extended$sIESs <- as.numeric(drop_out_df_extended$sIESs)

Fig2A <- ggscatter(drop_out_df_extended, x = "Total_IES_drop_out", y = "sIESs",
                   color = "#00468BFF", size = "Read_pairs_mapped", # Points color, shape and size
                   add = "reg.line",  # Add regression line
                   add.params = list(color = "#00AFBB", fill = "lightgray"), # Customize reg. line
                   conf.int = TRUE, # Add confidence interval
                   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                   cor.coeff.args = list(method = "pearson", label.x = 0.3, label.y = 450, label.sep = "\n"),
                   label = "sample", repel = TRUE,
                   xlim = c(0, 0.5),
                   ylim = c(0, 500),
                   xlab = "total IES dropout",
                   ylab = expression('somatic IESs'[''])) + 
  geom_point(data = drop_out_df_extended, 
             mapping = aes(x = Total_IES_drop_out, y = sIESs),
             color = "#00468BFF",
             size = 4 * drop_out_df_extended$Read_pairs_mapped / median(drop_out_df_extended$Read_pairs_mapped)) +
  theme_bw() + 
  theme(plot.background=element_rect(fill = "white"), 
        panel.background = element_rect(fill = 'white smoke'),
        axis.title.x = element_text(color = "black", size = 25, face = "plain", vjust=-2),
        axis.title.y = element_text(color = "black", size = 25, face = "plain", vjust= 0.5, angle = 90, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(color = "black", size = 20, face = "plain",  angle = 0),
        axis.text.y = element_text(color = "black", size = 20, face = "plain"),
        plot.margin = unit(c(1,2,1,1), "cm"))

Fig2A <- ggpar(Fig2A, 
               legend = "top",
               font.legend = c(15, "plain", "black"))
Fig2A

# Rearrange df
drop_out_df_extended$group <- "detected"
drop_out_df_extended1 <- drop_out_df_extended
drop_out_df_extended1$sample <- paste(drop_out_df_extended1$sample, rep("*", 12), sep = "")
drop_out_df_extended1$sIESs <- round(drop_out_df_extended$sIESs / (1-drop_out_df_extended$Total_IES_drop_out))
drop_out_df_extended1$group <- "inferred"
drop_out_df_extended_joined <- rbind(drop_out_df_extended, drop_out_df_extended1)
drop_out_df_extended_joined$group <- as.factor(drop_out_df_extended_joined$group)
drop_out_df_extended_joined$sample <- as.factor(drop_out_df_extended_joined$sample)

# Fig2B -------------------------------------------------------------------
# Build deviation plot
## Calculate the z-score of the sIES count
drop_out_df_extended_joined$IES_z[drop_out_df_extended_joined$group=="detected"] <- (drop_out_df_extended_joined$sIESs[drop_out_df_extended_joined$group=="detected"] - drop_out_df_extended_joined$sIESs[1])/sd(drop_out_df_extended_joined$sIESs[drop_out_df_extended_joined$group=="detected"])
drop_out_df_extended_joined$IES_z[drop_out_df_extended_joined$group=="inferred"] <- (drop_out_df_extended_joined$sIESs[drop_out_df_extended_joined$group=="inferred"] - drop_out_df_extended_joined$sIESs[1])/sd(drop_out_df_extended_joined$sIESs[drop_out_df_extended_joined$group=="inferred"])
drop_out_df_extended_joined$group <- factor(drop_out_df_extended_joined$group, levels = c("detected", "inferred"))

drop_out_df_extended_joined <- drop_out_df_extended_joined[drop_out_df_extended_joined$sample!="mcDNA*", ]
Fig2B <- ggbarplot(drop_out_df_extended_joined, x = "sample", y = "IES_z", 
                   fill = "group", 
                   color = "white",
                   palette = "jco",
                   sort.val = "asc",
                   label = c(sort(drop_out_df_extended_joined$sIESs[drop_out_df_extended_joined$group=="detected"]), 
                             sort(drop_out_df_extended_joined$sIESs[drop_out_df_extended_joined$group=="inferred"])),
                   lab.size = 3,
                   sort.by.groups = T, 
                   x.text.angle = 90,
                   ylab = "IES count (z-score)",
                   xlab = FALSE,
                   legend.title = "",
                   #label = round(drop_out_df_extended_joined$sIESs,1),
                   rotate = F,
                   ylim = c(-4, +4)) + 
  theme_bw() + 
  theme(plot.background=element_rect(fill = "white"), 
        panel.background = element_rect(fill = 'white smoke'),
        axis.title.x = element_text(color = "black", size = 25, face = "plain", vjust=-2),
        axis.title.y = element_text(color = "black", size = 25, face = "plain", vjust= 0.5, angle = 90, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(color = "black", size = 15, face = "plain",  angle = 90),
        axis.text.y = element_text(color = "black", size = 20, face = "plain"),
        plot.margin = unit(c(1,2,1,1), "cm"))

Fig2B <- ggpar(Fig2B,
               legend = "top",
               font.legend = c(15, "plain", "black"))

ggarrange(Fig2A, Fig2B, ncol = 2, nrow = 1, labels = c("A", "B"))

# lollipop chart ----------------------------------------------------------
ggdotchart(drop_out_df_extended_joined, x = "sample", y = "IES_z",
           color = "group",                              # Color by groups
           palette = c("#00AFBB", "#FC4E07"),            # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           group = "group",                                # Order by groups
           dot.size = 6,                                 # Large dot size
           label = round(drop_out_df_extended_joined$sIESs,1),  # Add IES counts as dot labels
           font.label = list(color = "white", size = 9, 
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr(),                       # ggplot2 theme
           ylim = c(-4, +4)) + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")
