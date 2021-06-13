# Description -------------------------------------------------------------
# Authors: Vitali et al. (2021)
# Table S4, basic experimental data
# IES accession numbers
# Empirical IRS estimates for bulk DNA-seq from mass culture
# Empirical IRS estimates for each scDNA-seq sample, at each time point  
# .........................................................................

# Libraries ---------------------------------------------------------------
library(tidyverse)

# Import data -------------------------------------------------------------
IES_tab <- read_tsv("./infiles/amitosis_IES_table.txt")
spec(IES_tab)

# Manipulate data ---------------------------------------------------------
## subset IES table and filter (aggressive)
## IES_COV > 20 in all 11 sc samples
## somatic IESs at D5 (IRS > 0.1)
## add ## mcDNA

## select columns of interest
IES_tab1 <- IES_tab %>%
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
            IRS_P25_sc4_D5 > 0.1) %>%
  print(n = Inf)

## Write to file
#write.table(IES_tab1, file = "../../ms/SciRep/Table_S4", 
#            quote = F, row.names = F)

# Assessment of IRS reliability after WGA ---------------------------------

### compute mean at D5 (IRS0)
IES_tab1 <- IES_tab1 %>%
  mutate(IRS_P25_sc0_D5_mean = rowMeans(across(ends_with("D5")))) %>%
  select(IES_ID, IRS_P25_sc0_D5_mean, everything(), -ends_with("D5"), -ends_with("D7"))

## reshape to long - Pivot1
IES_tab1_long <- pivot_longer(data = IES_tab1,
                              cols = matches("IRS"),
                              values_to = "IRS",
                              names_to = "sample",
                              names_prefix = "IRS_P25_")

## separate sample from day
IES_tab1_long1 <- separate(data = IES_tab1_long, 
                           col = "sample", 
                           into = c("sample", "day"), 
                           sep = 4)

## fix names
IES_tab1_long1$sample <- unlist(strsplit(IES_tab1_long1$sample, split = "_"))
IES_tab1_long1$day[grep("D5", IES_tab1_long1$day)] <- "D5"

## create combinatorial df
### create all possible D5-D10-D14 combinations
expanded <- expand_grid(IES_tab1_long1[IES_tab1_long1$day=="D5", ], 
                        IES_tab1_long1[IES_tab1_long1$day=="D10", ],
                        IES_tab1_long1[IES_tab1_long1$day=="D14", ],
                        .name_repair = "unique")

### keep combo for same IES
expanded$test <- ifelse(expanded$IES_ID...1 == expanded$IES_ID...5 & expanded$IES_ID...5 == expanded$IES_ID...9, TRUE, FALSE)
### rm duplicate colnames
expanded_D5_D10_D14 <- expanded[expanded$test==TRUE, -c(5, 9, 13)]

View(expanded_D5_D10_D14)

