# Description -------------------------------------------------------------
# Authors: Vitali et al. (2021)
# Simulate Somatic Assortment & Loss of Heterozygosity (H)
## Bioinformatic Simulation
# 1. Single locus simulation
# 2. Alternative Alleles encoded in binary digits
# 3. Ploidy level = 860n (Allen and Gibson 1972; Woodard et al. 1961)
# 4. Total number of copies conserved (constant ploidy)
# 5. Each daughter cell receives half copies
# 6. Selection-free environment
# .........................................................................

# Libraries ---------------------------------------------------------------
library(stringr)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(tictoc)

# Define Custom Functions --------------------------------------------------
## f1 amitosis = DNA Replication + Cell Division
amito <- function(m) {
  m1 <- paste(c(unlist(strsplit(m, split = "")), unlist(strsplit(m, split = ""))), collapse = "") # S phase
  D1 <- sample(unlist(strsplit(m1, split = "")), 860, replace = F) # Amitosis
  D2 <- c(rep(1, length(grep("1", unlist(strsplit(m1, split = "")))) - length(grep("1", D1))), 
          rep(0, length(grep("0", unlist(strsplit(m1, split = "")))) - length(grep("0", D1))))
  c(paste(D1, collapse = ""), paste(D2, collapse = ""))
}

## f2 convert to freq
conv2freq <- function(m) {
  f <- length(grep(1, unlist(strsplit(m, split = "")))) / nchar(m)
  as.vector(f)
}

# Simulate Daily Re-Isolation ---------------------------------------------
# 1. Single-cell bottlenecks
# 2. Growth Rate (GRate) = 4 div/24h
# 3. 2^10 replicate populations
# Total number of cells produced per day = 2^14 (2^10 * 2^4)
# .........................................................................
# Set number of replicates
REP = 2^10
# Set days 
DAY <- 50 # 50
# Growth Rate = 4 div/24h
div = seq(1, 4, 1)
# define final df
MAC_df_1 <- tibble(REP = 0,
                   DAY = 0,
                   GEN = 0, 
                   MAC = list(""))

for (r in 1:REP) {
  # define founder cell
  MAC0 <- paste(rep(c(1, 0), 430), collapse = "")
  # define storage df
  MAC_df <- tibble(REP = r, DAY = 1, GEN = 0, MAC = list(MAC0))
  # Reset frame row counter (w)
  w = 1
  # GEN counter (g)
  g = 0
  for (d in 1:DAY) { # isolation days
    for (i in seq_along(div)) { # divisions
      MAC_df[w + i, ] <- tibble(REP = r,
                                DAY = d,
                                GEN = g + i, 
                                MAC = list(lapply(unlist(MAC_df[(w - 1) + i, 4]), amito)))
    }
    # Pick a founder randomly from the 2^4 cells (single-cell bottleneck)
    MAC0 <- sample(unlist(MAC_df[w + i, 4]), 1, replace = F)
    MAC_df[w + i + 1, ] <- tibble(REP = r, DAY = d + 1, GEN = MAC_df$GEN[w + i], MAC = list(MAC0))
    # increment row counter
    w = w + 5
    # increment GEN counter
    g = g + 4
  }
  # Bind REPs
  MAC_df_1 <- bind_rows(MAC_df_1, MAC_df)
}

## Remove first entry (0, 0, 0)
MAC_df_1 <- MAC_df_1[-1, ]

## Select generations with 16 cells (daily offspring)
MAC_df_1 <- MAC_df_1 %>%
  mutate(cells = unlist(map(MAC_df_1$MAC, ~length(unlist(.))))) %>% 
  filter(cells == 16)

#MAC_df_1 %>%
  #print(n = Inf)

## Reshape to long (separate cells into single entries)
MAC_df_long <- tibble(REP = 1, DAY = 0, GEN = 0, MAC = unlist(MAC_df_1[1, 4]))
for (i in 2:nrow(MAC_df_1)) {
  MAC_df_long <- bind_rows(MAC_df_long,
                           tibble(MAC_df_1[i, 1],
                                  MAC_df_1[i, 2],
                                  MAC_df_1[i, 3],
                                  MAC = unlist(unlist(MAC_df_1[i, 4]))))
}

## Convert to Frequency
MAC_df_long$IRS <- as.vector(unlist(lapply(MAC_df_long$MAC, conv2freq)))
# Check distribution for last GEN
summary(MAC_df_long$IRS[MAC_df_long$GEN == DAY * length(div)])

# Visualize Allele Freq Distribution DR------------------------------------
## e.g. IES+ alleles = IRS distribution
## GEN subset
gen_sub <- seq(4, DAY * length(div), 60) # 16; 8
MAC_df_long_subset <- MAC_df_long %>%
  filter(GEN %in% gen_sub)

# GEN to factor
MAC_df_long_subset$GEN <- as.factor(MAC_df_long_subset$GEN)
# set palette
pal = pal_lancet("lanonc")(length(gen_sub))

ggplot(MAC_df_long_subset, aes(x = IRS, color = GEN, ..density..)) +
  geom_density() +
  xlim(c(0, 1.0)) +
  scale_color_manual(values = pal) +
  labs(title="", x = expression("IES"^'+'*"copies (fraction)"), 
       y = "kernel density") +
  theme_bw() + 
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = 'white smoke'),
        axis.title.x = element_text(color = "black", size = 25, face = "plain", vjust=-1),
        axis.title.y = element_text(color = "black", size = 25, face = "plain", vjust= 0.5, 
                                    angle = 90, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(color = "black", size = 15, face = "plain"),
        axis.text.y = element_text(color = "black", size = 15, face = "plain"),
        plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "top")

# Rate of Somatic Assortment ----------------------------------------------
## Change of IRS SD across generations across all runs
SD_change <- MAC_df_long %>%
  filter(GEN != 0) %>% 
  group_by(GEN) %>%
  summarise(SD = sd(IRS)) %>%
  print(n = Inf)
### Add group label
SD_change$Simulation <- "Daily re-isolation (n=1024)"

# Define x and y
y <- SD_change$SD
x <- SD_change$GEN

# Model-1, Quadratic model
m1 <- lm(y ~ poly(x, 2))
summary(m1)
confint(m1)
cor(y,predict(m1)) # goodness of fit estimate

# plot SD values across GEN
plot(x, y,  pch = 16, 
     col = "dark red", xlab = "GEN",
     ylab = expression('SD'[(IRS)]),
     main = "IRS0 = 0.5")

# plot  fitted line, Model-1
lines(x, predict(m1, data.frame(x = x)), lty = 2, col = "blue", lwd = 2)

# Predict values and C.I.99
predicted.intervals <- predict(m1, data.frame(x = x), interval = 'confidence',
                               level = 0.95)


# Validate Mathematical Modeling -------------------------------------------
## Figure S5
## Change of IRS SD across generations for small samples chosen at random (n=4)
SD_change1 <- MAC_df_long[MAC_df_long$REP == sample(unique(MAC_df_long$REP), 4), ] %>%
  filter(GEN != 0) %>% 
  group_by(GEN) %>%
  summarise(SD = sd(IRS)) %>%
  print(n = Inf)
### Add group label
SD_change1$Simulation <- "Daily re-isolation (n=4)"

## Load SENES core module function
load("./fun/senes_core.RData")

## set up frames
hetero_df <- list()
hetero_df_all <- list()

## Define parameters
k = 860 # number of elements drawn (860)
Chr = 1 # number of somatic chromosomes
N = 2*k*Chr # total number of elements to draw from, number of elements after phase S (2x860)
GEN <- seq(2, 200, 1) # set number of generations (GEN = 200, full clonal cycle)
IRS0 <- seq(0.5, 0.5, 0.1) # Vector of IRS0
spread = "" # standard deviation

## execute core module
tic("simulation time")
hetero_df_all <- senes_core(k, Chr, N, IRS0, GEN, hetero_df, hetero_df_all)
toc()

## convert to long format (P)
math_sim_long <- as.data.frame(matrix(ncol = ncol(hetero_df_all)+1, nrow = 0))
names(math_sim_long) <- c(names(hetero_df_all), "Alleles_n")
for (i in 1:length(hetero_df_all$GEN)) {
  temp <- cbind.data.frame(hetero_df_all[i, 1:5], unlist(hetero_df_all$P_dist[i]))
  temp$Alleles_n <- seq(1, (k+1), 1)-1 # add x
  names(temp) <- names(math_sim_long)
  math_sim_long <- rbind(math_sim_long, temp)
}

## Manipulate frame
### to factor
math_sim_long$GEN <- as.factor(math_sim_long$GEN)
math_sim_long$IRS0 <- as.factor(math_sim_long$IRS0)
### 
SD_change_math <- math_sim_long %>%
  filter(GEN != 0, !duplicated(GEN)) %>% 
  select(GEN, sd)
### convert GEN to numeric
SD_change_math$GEN <- as.numeric(as.character(SD_change_math$GEN))
### rename variable sd -> SD
colnames(SD_change_math)[2] <- "SD"
### Add group label
SD_change_math$Simulation <- "Mathematical"

## combine frames
SD_change_combo <- rbind(SD_change, SD_change1, SD_change_math)

## Compare simulations
### set palette
pal = pal_lancet("lanonc")(length(gen_sub))
### plot
ggplot() +
  geom_point(SD_change_combo[SD_change_combo$Simulation == "Daily re-isolation (n=1024)", ],
             mapping = aes(GEN, SD, color = Simulation)) +
  geom_line(SD_change_combo[SD_change_combo$Simulation == "Mathematical", ],
            mapping = aes(GEN, SD, color = Simulation), 
            size = 0.8, linetype = 2) +
  geom_line(SD_change_combo[SD_change_combo$Simulation == "Daily re-isolation (n=4)", ],
            mapping = aes(GEN, SD, color = Simulation)) +
  scale_color_manual(values = pal) +
  ylab(expression('SD'[(IRS)])) +
  xlab("gen (asexual)") +
  theme_bw(base_size = 20) +
  theme(legend.position="top", 
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))

# Rate of loss of Heterozygosity ------------------------------------------
# Heterozygosity (H), fraction of nuclei with 0 < IRS < 1
H_change <- MAC_df_long %>% 
  mutate(H_lgc = ifelse(!near(IRS, 0) | !near(IRS, 1), TRUE, FALSE)) %>% 
  group_by(GEN) %>%
  summarise(H = mean(H_lgc)) %>%
  print(n = Inf)

ggplot(H_change, aes(GEN, H)) +
  geom_line() +
  theme_bw()

# fit nls to data
x1 = H_change$GEN
y1 = H_change$H

# Model-1, Cubic model
M1 <- lm(y1 ~ poly(x1, 3))
summary(M1)
confint(M1)
cor(y1, predict(M1)) # goodness of fit estimate

# plot H values
plot(x1, y1,  pch = 16, ylim = c(0, 1), col = "gray")
# plot  fitted line, Model-1
lines(x, predict(m1, data.frame(x = x1)), lty = 2, col = "dark red", lwd = 2)

# Predict values and C.I.99
predicted.intervals <- predict(M1, data.frame(x = x1), interval = 'confidence',
                               level = 0.95)

lines(H_change$GEN, predicted.intervals[, 1], col = 'red', lwd = 2)
lines(H_change$GEN, predicted.intervals[, 2], col = 'black', lwd = 1)
lines(H_change$GEN, predicted.intervals[, 3], col = 'black', lwd = 1)

# Simulate Mass Culture ---------------------------------------------------
# 1. 4 div/24h
# 2. bottlenecks of 1024 cells
# .........................................................................
## Create df for 12 Generations (2^12 = 4096 cells)
## Expand culture from single cell
# Set days 
DAY <- 50
# Growth Rate = 4 div/24h
div = seq(1, 4, 1)
# divisions to seed first mass culture
div_seed = seq(1, 12, 1)
# Bottleneck (BN)
BN = 2^10
# define founder cell
MAC0 <- paste(rep(c(1, 0), 430), collapse = "")
# define storage df
MAC_df_ms <- tibble(GEN = 0, MAC = list(MAC0))

# Simulate seed culture
for (i in seq_along(div_seed)) { # divisions
  MAC_df_ms[i + 1, ] <- tibble(GEN = i,
                               MAC = list(lapply(unlist(MAC_df_ms[i, grep("MAC", colnames(MAC_df_ms))]), amito)))
}

# Add DAY
MAC_df_ms$DAY <- c(0, rep(1, 4), rep(2, 4), rep(3, 4))
# Move DAY to first column
MAC_df_ms <- MAC_df_ms %>% 
  relocate(DAY, before = GEN)
# Collect inocolum (1024 cells)
MAC0 <- sample(unlist(MAC_df_ms[nrow(MAC_df_ms), grep("MAC", colnames(MAC_df_ms))]), BN, replace = FALSE)
# Bind founders to frame
MAC_df_ms <- bind_rows(MAC_df_ms, tibble(DAY = 4, GEN = 12, MAC = list(MAC0)))
# Set frame row counter (w)
w = nrow(MAC_df_ms)
# GEN counter (g)
g = as.numeric(MAC_df_ms[nrow(MAC_df_ms), grep("GEN", colnames(MAC_df_ms))])

for (d in seq_along(4:DAY)) { # mass culture days
  for (i in seq_along(div)) { # divisions per day
    MAC_df_ms[w + i, ] <- tibble(DAY = seq(4, DAY)[d],
                              GEN = g + i, 
                              MAC = list(lapply(unlist(MAC_df_ms[(w - 1) + i, grep("MAC", colnames(MAC_df_ms))]), amito)))
  }
  # Pick founders randomly to seed next culture (1024-cell bottleneck)
  MAC0 <- sample(unlist(MAC_df_ms[w + i, grep("MAC", colnames(MAC_df_ms))]), BN, replace = FALSE)
  MAC_df_ms[w + i + 1, ] <- tibble(DAY = seq(4, DAY)[d] + 1, GEN = MAC_df_ms$GEN[w + i], MAC = list(MAC0))
  # increment row counter
  w = w + 5
  # increment GEN counter
  g = g + 4
}

## Select generations with 2^14 cells (daily offspring from 2^10 inocolum)
MAC_df_ms <- MAC_df_ms %>%
  mutate(cells = unlist(map(MAC_df_ms$MAC, ~length(unlist(.))))) %>% 
  filter(cells == 2^14)

#MAC_df_ms %>%
  #print(n = Inf)

## Reshape to long (separate cells into single entries)
MAC_df_ms_long <- tibble(DAY = 4, GEN = 16, MAC = unlist(MAC_df_ms[1, grep("MAC", colnames(MAC_df_ms))]))
for (i in 2:nrow(MAC_df_ms)) {
  MAC_df_ms_long <- bind_rows(MAC_df_ms_long,
                           tibble(MAC_df_ms[i, 1],
                                  MAC_df_ms[i, 2],
                                  MAC = unlist(unlist(MAC_df_ms[i, grep("MAC", colnames(MAC_df_ms))]))))
}

## Convert to Frequency
MAC_df_ms_long$IRS <- as.vector(unlist(lapply(MAC_df_ms_long$MAC, conv2freq)))
# Check distribution for last GEN
summary(MAC_df_ms_long$IRS[MAC_df_ms_long$GEN == DAY * length(div)])

# Visualize Allele Freq Distribution MC------------------------------------
## e.g. IES+ alleles = IRS distribution
## GEN subset
gen_sub <- seq(16, DAY * length(div), 36)
MAC_df_ms_long_subset<- MAC_df_ms_long %>%
  filter(GEN %in% gen_sub)

# GEN to factor
MAC_df_ms_long_subset$GEN <- as.factor(MAC_df_ms_long_subset$GEN)
# set palette
pal = pal_lancet("lanonc")(length(gen_sub))

ggplot(MAC_df_ms_long_subset, aes(x = IRS, color = GEN, ..density..)) +
  geom_density() +
  xlim(c(0, 1.0)) +
  scale_color_manual(values = pal) +
  labs(title="", x = expression("IES"^'+'*"copies (fraction)"), 
       y = "kernel density") +
  theme_bw() + 
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = 'white smoke'),
        axis.title.x = element_text(color = "black", size = 25, face = "plain", vjust=-1),
        axis.title.y = element_text(color = "black", size = 25, face = "plain", vjust= 0.5, 
                                    angle = 90, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(color = "black", size = 15, face = "plain"),
        axis.text.y = element_text(color = "black", size = 15, face = "plain"),
        plot.margin = unit(c(1,1,1,1), "cm"), legend.position = "none")

# Change of IRS SD across generations -------------------------------------
SD_change_mc <- MAC_df_ms_long %>%
  group_by(GEN) %>%
  summarise(SD = sd(IRS)) %>%
  print(n = Inf)

ggplot(SD_change_mc[-1, ], aes(GEN, SD)) +
  geom_line() +
  theme_bw()
