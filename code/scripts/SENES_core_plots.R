# Description -------------------------------------------------------------
# Authors: Vitali et al. (2021) 
# Simulate Somatic Assortment & Loss of Heterozygosity (H)
## Mathematical Simulation
# 1. Single locus simulation
# 2. Hypergeometric distribution to model first generation
# 3. Ploidy level = 860n (Allen and Gibson 1972; Woodard et al. 1961)
# 4. Total number of copies conserved (constant ploidy)
# 5. Each daughter cell receives half copies
# 6. Selection-free environment
# .........................................................................

# Libraries ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggsci)
library(tictoc)

# SENES define core module function ---------------------------------------
senes_core <- function(k, Chr, N, IRS0, GEN, hetero_df, hetero_df_all) {
  for (w in 1:length(IRS0)) {
    m = round(2*IRS0[w]*k) # number of success states, IES+ copies available (integer)
    input = m / 2 / k # input_ratio
    if (input > 0.5) {
      m = round(2*(1-IRS0[w])*k)
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
      #y <- Reduce("+", d) # neat but a tad slower than the for loop
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
      # ................................
      spread <-  sqrt(sum(x1^2 * y2) - sum(x1 * y2)^2) # calculate sd for num of IES+ copies
      hetero_df[[j]] <- cbind(j+1, IRS0[w], y2[1]+y2[k+1], sum(y2[2:(k)]), spread[[1]] / k, list(list(y2)))
      print(c("IRS0:", IRS0[w], "|", "GEN:", GEN[j]), quote = F)
      
    }
    
    hetero_df_all[[w]] <- hetero_df
    
  }
  
  hetero_df_all <- lapply(hetero_df_all, function(x) do.call(rbind, x))
  hetero_df_all <- do.call(rbind.data.frame, hetero_df_all)
  names(hetero_df_all) <- c("GEN","IRS0", "Homoz", "H", "sd", "P_dist")
  hetero_df_all <- unnest(hetero_df_all, cols = 1:5)
  return(hetero_df_all)
  
}

# set up dataframes
hetero_df <- list()
hetero_df_all <- list()

# Define parameters -------------------------------------------------------
k = 860 # number of elements drawn (860)
Chr = 1 # number of somatic chromosomes
N = 2*k*Chr # total number of elements to draw from, number of elements after phase S (2x860)
GEN <- seq(2, 200, 1) # set number of generations, e.g. D9 (GEN = 200, full clonal cycle)
IRS0 <- seq(0.1, 0.9, 0.1) # Vector of IRS0 to calculate the rate of SA and H loss
#IRS0 <- 0.5
spread = "" # standard deviation

# Execute core module function --------------------------------------------
tic("simulation time")
hetero_df_all <- senes_core(k, Chr, N, IRS0, GEN, hetero_df, hetero_df_all)
toc()

# Manipulate data frame ---------------------------------------------------
## convert to long format (P) and combine
reduced_all <- as.data.frame(matrix(ncol = ncol(hetero_df_all)+1, nrow = 0))
names(reduced_all) <- c(names(hetero_df_all), "Alleles_n")
for (i in 1:length(hetero_df_all$GEN)) {
  temp <- cbind.data.frame(hetero_df_all[i, 1:5], unlist(hetero_df_all$P_dist[i]))
  temp$Alleles_n <- seq(1, (k+1), 1)-1 # add x
  names(temp) <- names(reduced_all)
  reduced_all <- rbind(reduced_all, temp)
}

## to factor
reduced_all$GEN <- as.factor(reduced_all$GEN)
reduced_all$IRS0 <- as.factor(reduced_all$IRS0)

## write long df to file
#write.table(reduced_all, paste("./senes", "chr", Chr, "k", k, "gen", 
#                               max(GEN), sep = "_"), sep = "\t" , quote = FALSE)


# Plot - p(X) IRS0 = 0.1, 200 GEN -----------------------------------------
reduced_all_sub <- reduced_all %>%
  filter(GEN %in% c(2, 10, 20, 50, 200), IRS0 == 0.1)

## set palette
pal = pal_lancet("lanonc")(5)

Fig1D <- ggplot(reduced_all_sub, aes(x = Alleles_n / k, P_dist, color = GEN)) +
  geom_line() +
  xlim(c(0, 0.3)) +
  scale_color_manual(values = pal) +
  labs(title="", x = expression("IES"^'+'*"copies (fraction)"), 
       y = expression("P"['(IES'^'+'*', gen)'])) +
  theme_bw() + 
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = 'white smoke'),
        axis.title.x = element_text(color = "black", size = 25, face = "plain", vjust=-1),
        axis.title.y = element_text(color = "black", size = 25, face = "plain", vjust= 0.5, 
                                    angle = 90, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(color = "black", size = 15, face = "plain"),
        axis.text.y = element_text(color = "black", size = 15, face = "plain"),
        plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "none")

ggpar(Fig1D, legend = "top", 
      font.legend = c(15, "plain", "black"))


# plot - p(X) distribution ------------------------------------------------
# Facet each plot by IRS
## reduce df to a subset of GENs (filter)
GEN_sub <- c(seq(GEN[1], GEN[length(GEN)], GEN[length(GEN)]/4),  GEN[length(GEN)])
reduced_all <- as_tibble(reduced_all) %>%
  filter(GEN %in% GEN_sub)
## set color palette
pal = colorRampPalette(pal_jama()(7))(length(levels(reduced_all$GEN)))
## set y limit to improve visualization (for coord_cartesian)
roof = reduced_all %>% 
  as_tibble() %>%
  filter(Alleles_n > 0 & Alleles_n < k)
roof = max(roof$P_dist) * 1.2

ggplot(reduced_all, aes(x = Alleles_n, P_dist, color = GEN)) +
  geom_line() +
  xlim(c(0, k)) +
  scale_color_manual(values = pal) +
  coord_cartesian(ylim=c(0, roof)) +
  facet_wrap(~ IRS0, ncol = 3) + # facet_grid(GEN ~ IRS0) +
  theme_bw() + 
  theme(plot.background=element_rect(fill = "white"), 
        panel.background = element_rect(fill = 'white smoke'),
        axis.title.x = element_text(color = "black", size = 25, face = "plain", vjust=-2),
        axis.title.y = element_text(color = "black", size = 25, face = "plain", vjust= 0.5, 
                                    angle = 90, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(color = "black", size = 10, face = "plain"),
        axis.text.y = element_text(color = "black", size = 10, face = "plain"),
        plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "none")


# Plot - Rate of Somatic Assortment ---------------------------------------
hetero_df_all$IRS0 <- as.factor(hetero_df_all$IRS0)
Fig1B <- ggline(hetero_df_all, x = "GEN", y = "sd", 
                color = "IRS0", 
                palette = "lancet",
                plot_type = "l",
                size = 0.5,
                numeric.x.axis = T,
                ylim = c(0, max(hetero_df_all$sd)*1.5), # c(0, 0.2)
                xlim = c(1, length(GEN)),
                xlab = "gen (asexual)",
                ylab = expression('SD'['(IRS)'])) +
  #title = "Rate of Somatic Assortment") + 
  theme_bw() + 
  theme(plot.background=element_rect(fill = "white"), 
        panel.background = element_rect(fill = 'white smoke'),
        axis.title.x = element_text(color = "black", size = 25, face = "plain", vjust=-2),
        axis.title.y = element_text(color = "black", size = 25, face = "plain", vjust= 0.5, 
                                    angle = 90, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(color = "black", size = 15, face = "plain"),
        axis.text.y = element_text(color = "black", size = 15, face = "plain"),
        plot.margin = unit(c(1,1,1,1), "cm"))

Fig1B <- ggpar(Fig1B, legend = "top", 
               font.legend = c(15, "plain", "black"), 
               xticks.by = round(length(GEN)/10))


# Plot - Loss of Heterozygosity -------------------------------------------
## heterozygosity, fraction of nuclei with 0 < IRS < 1
hetero_df_all$IRS0 <- as.factor(hetero_df_all$IRS0)
Fig1C <- ggline(hetero_df_all, x = "GEN", y = "H", 
                color = "IRS0", 
                palette = "lancet",
                plot_type = "l",
                size = 0.5,
                numeric.x.axis = T,
                ylim = c(0, 1),
                xlim = c(1, GEN_sub[length(GEN_sub)]),
                xlab = "gen (asexual)",
                ylab = expression('H'['(IES)'])) +
  #title = "Loss of heterozygosity (H)") + 
  theme_bw() + 
  theme(plot.background=element_rect(fill = "white"), 
        panel.background = element_rect(fill = 'white smoke'),
        axis.title.x = element_text(color = "black", size = 25, face = "plain", vjust=-2),
        axis.title.y = element_text(color = "black", size = 25, face = "plain", vjust= 0.5, angle = 90, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(color = "black", size = 12, face = "plain"),
        axis.text.y = element_text(color = "black", size = 15, face = "plain"),
        plot.margin = unit(c(1,1,1,1), "cm"))

Fig1C <- ggpar(Fig1C, legend = "top", 
               font.legend = c(15, "plain", "black"), 
               xticks.by = round(length(GEN)/4))


# Compare SENES with Preer Genet. Res. (1976) -----------------------------
# haploid model
### read in pre-calculated data (SENES)
senes_chr_1_k_860_gen_250 <- read.table("./infiles/senes_chr_1_k_860_gen_250", h = T)
senes_chr_1_k_860_sub <- senes_chr_1_k_860_gen_250 %>% 
  filter(GEN %in% seq(50, 250, 50), Alleles_n == 430)

Preer_Tab5 <- c(0.169, 0.238, 0.289, 0.331, 0.368)
data.frame(GEN = senes_chr_1_k_860_sub$GEN, 
           senes_sd = round(senes_chr_1_k_860_sub$sd*2, digits = 3), 
           Preer1976_sd = Preer_Tab5)

### fraction of cells with IRS >= 0.85 after 200 generations
senes_chr_1_k_860_gen_250 %>% 
  filter(GEN == 200, Alleles_n >= 0.85 * 860) %>% 
  summarise(Cumul = sum(P_dist))

## chromosomal model
### read in pre-calculated data SENES
senes_chr_43_k_860_gen_250 <- read.table("./infiles/senes_chr_43_k_860_gen_250", h = T)
senes_chr_43_k_860_sub <- senes_chr_43_k_860_gen_250 %>% 
  filter(GEN %in% seq(50, 250, 50), Alleles_n == 430)

Preer_Tab3 <- c(0.240, 0.339, 0.415, 0.479, 0.536)
data.frame(GEN = senes_chr_43_k_860_sub$GEN, 
           senes_sd = round(senes_chr_43_k_860_sub$sd*2, digits = 3), 
           Preer1976_sd = Preer_Tab3)

### fraction of cells with IRS >= 0.85 after 200 generations
senes_chr_43_k_860_gen_250 %>% 
  filter(GEN == 200, Alleles_n >= 0.85 * 860) %>% 
  summarise(Cumul = sum(P_dist))
