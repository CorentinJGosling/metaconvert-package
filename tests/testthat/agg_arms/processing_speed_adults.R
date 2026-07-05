library(readxl)
library(dplyr)
library(tidyr)
library(tidyverse)

# data <- read_excel("C:/Users/cc1d24/OneDrive - University of Southampton/UoS team/Projects/Robert/Analysis Robert/28_01_26_Database_Neurocog_Onlyrelevant.xlsx")


# data <- read_excel("C:/Users/coren/Documents/sideprojet/christos-NMA/28_01_26_Database_Neurocog_Onlyrelevant.xlsx")
data <- read_excel("28_01_26_Database_Neurocog_Onlyrelevant.xlsx")

names(data) <- data[2, ]
data <- data[-c(1,2), c(1:59, 61:66)]
names(data) <- make.unique(names(data))

data <- data %>% filter(!(cog_domain %in% c("DO NOT USE", "unclear", "NOT AVAILABLE", "`"))) 


########################### Child/adolescents ##################################


########### Reaction time/processing speed outcome ##########

speed <- data %>% filter(cog_domain == "Processing speed") %>%
  filter(`Pediatric_or_adult_sample?` == 4)
speed <- speed[, c(3, 6, 12:18, 28:65)]




###################### Endpoind ####################

endp <- speed %>%
  select(7, 3, 9, 10, 12, 21, 22, 23, 34:36, 32) %>%
  mutate(
    Neurocog_bas_nitt_12 = ifelse(Neurocog_end_mean_12 != "$", Neurocog_bas_nitt_12, "$"),
    Neurocog_end_mean_12 = ifelse(Neurocog_end_mean_12 != "$", Neurocog_end_mean_12, Neurocog_end_mean_compl_12),
    Neurocog_end_sd_12 = ifelse(Neurocog_end_sd_12 != "$", Neurocog_end_sd_12,   Neurocog_end_sd_compl_12),
    Neurocog_bas_nitt_12 = ifelse(Neurocog_bas_nitt_12 != "$", Neurocog_bas_nitt_12, Neurocog_N_compl_12),
    Neurocog_bas_nitt_12 = as.numeric(Neurocog_bas_nitt_12), Neurocog_end_sd_12 = as.numeric(Neurocog_end_sd_12),
    Neurocog_end_mean_12 = as.numeric(Neurocog_end_mean_12)) %>%
  select("Study ID", Task, measure, Compound_list, better_if_small_or_large, Neurocog_bas_nitt_12,                   
         Neurocog_end_mean_12, Neurocog_end_sd_12) %>%
  filter(!is.na(Neurocog_end_mean_12))


# Impute some sd's from se's
endp[endp$Neurocog_bas_nitt_12 == 43, 8] <- sqrt(43)*0.4
endp[endp$Neurocog_bas_nitt_12 == 51, 8] <- sqrt(51)*0.4


# Exclude Taylor 2001 (one of two for not having sd)
endp <- endp %>% filter(!is.na(Neurocog_end_sd_12))



# Make it wide 
endp_wide <- endp %>%
  full_join(endp, 
            by = c("Study ID", "Task", "measure", "better_if_small_or_large"),
            suffix = c("_arm1", "_arm2")) %>%
  filter(Compound_list_arm1 < Compound_list_arm2) %>%
  select(`Study ID`, Task, measure, better_if_small_or_large,
         Compound_list_arm1, Compound_list_arm2, Neurocog_bas_nitt_12_arm1, Neurocog_bas_nitt_12_arm2, 
         Neurocog_end_mean_12_arm1, Neurocog_end_mean_12_arm2, everything()) %>% 
  rename(sample1 = Neurocog_bas_nitt_12_arm1, sample2 = Neurocog_bas_nitt_12_arm2, mean1 =  Neurocog_end_mean_12_arm1, 
         mean2 =  Neurocog_end_mean_12_arm2, sd1 = Neurocog_end_sd_12_arm1, sd2 = Neurocog_end_sd_12_arm2,
         treat1 = Compound_list_arm1, treat2 = Compound_list_arm2)


rm(endp)



# Calculate the effect sizes
library(metafor)
es.endp <-  escalc(measure = "SMD", n1i = sample1, n2i = sample2, 
                   m1i = mean1, m2i = mean2, sd1i = sd1, sd2i = sd2, 
                   data = endp_wide) %>% 
  select(1, 2, 3, 4, 5, 6, 13, 14)


# Now change the association in the lower better
es.endp <- es.endp %>%
  mutate(yi = ifelse(better_if_small_or_large == "lower", yi * -1, yi)) %>% select(-better_if_small_or_large)

rm(endp_wide)





## !!!!!!!Corentin here is the problem that I told you about (Taylor trials)!!!!!

library(metaConvert)
metacon <- aggregate_df(es.endp, dependence = "outcomes", cor_unit = 0.5, 
                        agg_fact = "Study.ID", es = "yi", se = "vi")

es.end <- es.endp[, c(1, 4, 5)] %>% distinct() %>% left_join(metacon[, 2:4], by = "Study.ID")
rm(es.endp, metacon)

