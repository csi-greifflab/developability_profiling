# setup ####
library(tidyverse)
library(data.table)
library(ggpubr)
library(gghalves)
library(moments)

AbStruc_native_developability_final <- fread("data/AbStruc_native_developability_final.csv")[, -1] # WT structure DP
mutated_struc_dp <- read_csv("data/AbStruc_mutated_cdrs_developability_igfold_63K.csv")
wt_struc_dp <- readRDS("data/igfold_sDPs_WT_all.rds")
mutated_seq_dp <- fread("data/AbChain_mutated_full_developability.csv")

# order and clean datasets ####
filter_cols <- function(df) {
  df %>%
    select(-c(
      ends_with(c("gene", "percentextcoef", 'content')),
      starts_with(c("aaSeq", "fr", "cdr")),
      starts_with(c(
        "database",
        "chain",
        "AbChain_length",
        "AbChain_cysbridges_percenextcoef",
        "AbChain_1_charge", "AbChain_2_charge", "AbChain_3_charge",
        "AbChain_5_charge", "AbChain_6_charge", "AbChain_8_charge",
        "AbChain_9_charge", "AbChain_10_charge", "AbChain_11_charge",
        "AbChain_13_charge", "AbChain_14_charge",
        'AbStruc_steric_clashes' #Didn't work
      ))
    )) %>%
    na.omit() %>% 
    mutate(rowid = as.character(rowid))
}

mutated_seq_dp <- mutated_seq_dp %>% 
  as_tibble() %>%
  filter_cols()

mutated_struc_dp <- mutated_struc_dp %>%
  as_tibble() %>% 
  filter_cols() %>%
  filter(rowid != '2099484') # 2099484 no succesful igFold DP predictions

# Get WT rowid
struc_rowids <- mutated_struc_dp %>%
  pull(rowid) %>%
  unique()

# Extract metadata and add to mutated dataset
wt_meta <- mutated_seq_dp %>%
  filter(mutation_type == "none") %>%
  select(c(rowid, isotype, species)) %>% 
  mutate(rowid = as.character(rowid))

# Extract DP of wt and filter unwanted columns
wt_strucs <- wt_meta %>%
  filter(rowid %in% struc_rowids) %>%
  mutate(mutation_type = "none") %>%
  left_join(wt_struc_dp) %>%
  filter_cols()
  

regions = mutated_seq_dp %>%
  select(rowid, mutation_type, mutation_region)

mutated_struc_dp <- mutated_struc_dp %>%
  left_join(wt_meta) %>%
  bind_rows(wt_strucs) %>% # add WT data
  left_join(regions) %>% # add region data
  filter(mutation_region %in% c('CDR1', 'CDR2', 'CDR3', 'none')) # Keep only CDR mutants and wildtype

# set isotype and mutation_region as factors and manually define levels
mutated_seq_dp = mutated_seq_dp %>%
  mutate(mutation_region = as.factor(mutation_region)) %>%
  mutate(isotype = as.factor(isotype)) %>%
  mutate(species = as.factor(species)) %>%
  mutate(mutation_region = fct_relevel(mutation_region, c("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"))) %>%
  mutate(isotype = fct_relevel(isotype, c("IgD", "IgM", "IgG", "IgA", "IgE"))) %>%
  filter(species == 'human') %>%
  filter(!isotype %in% c('IgL', 'IgK'))

mutated_struc_dp = mutated_struc_dp %>%
  left_join(mutated_seq_dp) %>%
  mutate(mutation_region = as.factor(mutation_region)) %>%
  mutate(isotype = as.factor(isotype)) %>%
  mutate(species = as.factor(species)) %>%
  mutate(mutation_region = fct_relevel(
    mutation_region,
    c("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")
  )) %>%
  mutate(isotype = fct_relevel(isotype, c(
    "IgD", "IgM", "IgG", "IgA", "IgE"
  ))) %>%
  filter(species == 'human') %>%
  filter(!isotype %in% c('IgL', 'IgK'))

rm(wt_meta, wt_strucs, struc_rowids, AbStruc_native_developability_final)

save(mutated_seq_dp, mutated_struc_dp, file = '/storage/jahnz/R/developability/data/sensitivity/human_heavy_mutations_filtered.RData')
load('/storage/jahnz/R/developability/data/sensitivity/human_heavy_mutations_filtered.RData')

# calculate distances ####
# mutated_seq_dp = mutated_seq_dp %>%
#   filter(mutation_region %in% c('CDR1','CDR2','CDR3'))
mutated_struc_dp = mutated_struc_dp %>% select(-c(starts_with('AbChain')))

scale_manual_fun <- function(x) {
  (x - mean(x)) / sd(x)
}

scale_data = function(df) {
  df %>%
    mutate_at(scale_manual_fun, .vars = vars(starts_with(c("Ab"))))
}

# Get DP of WT strucs
extract_original <- function(df) {
  df %>%
    filter(mutation_type == "none") %>%
    transmute(rowid = rowid, across(starts_with("Ab"), .names = "wt_{col}"))
}

# gets difference between mutated and wt DP value
get_distance <- function(df) {
  df %>%
    filter("mutation_type" != "none") %>%
    left_join(extract_original(df)) %>% 
    scale_data() %>%
    mutate(abs(across(starts_with("Ab"), .names = "{col}_dist") - across(starts_with("wt")))) %>%
    select(c("rowid", starts_with('mutation'), "isotype", "species", ends_with("dist")))
}

rmsd = read_csv('/storage/jahnz/R/developability/data/sensitivity/sensitivity_RMSD_WT_mutants_IgFold.csv')
rmsd = rmsd %>%
  mutate(rmsd = as.numeric(rmsd)) %>%
  mutate(rowid = as.character(ref)) %>%
  select(c(rowid, mobile, rmsd)) %>%
  separate(mobile, c(NA,'mutation_type'), sep='_')

mutated_seq_dp_dist = mutated_seq_dp %>% 
  get_distance() 

mutated_struc_dp_dist = mutated_struc_dp %>% 
  get_distance() %>%
  left_join(rmsd)

save(mutated_seq_dp_dist,
     mutated_struc_dp_dist,
     file = '/storage/jahnz/R/developability/Figures/Thesis/human_heavy_mutated_dp_distances.RData')



# Calculate sensitivity metrics ####

# load dp distances
load('/storage/jahnz/R/developability/Figures/Thesis/human_heavy_mutated_dp_distances.RData')


# Overall sensitivity
overall_seq_kurtosis = mutated_seq_dp %>%
  filter(mutation_region != 'none') %>%
  scale_data() %>%
  pivot_longer(-c(-starts_with('Ab')),
               names_to = 'parameter') %>%
  group_by(rowid, parameter) %>%
  summarise(kurtosis = kurtosis(value) - 3,
            range = diff(range(value))) %>%
  ungroup() %>%
  select(rowid, parameter, kurtosis, range)


overall_seq_sensitivity = mutated_seq_dp_dist %>%
  filter(mutation_region != 'none') %>%
  pivot_longer(-c(rowid, mutation_type, mutation_region, isotype, species),
               values_to = "distance",
               names_to = "parameter") %>%
  group_by(rowid, parameter) %>%
  summarise(ASD = median(distance)) %>%
  ungroup() %>%
  mutate(parameter = gsub("_dist","",parameter)) %>% 
  left_join(overall_seq_kurtosis)

overall_struq_kurtosis = mutated_struc_dp %>%
  filter(mutation_region != 'none') %>%
  scale_data() %>%
  pivot_longer(-c(-starts_with('Ab')),
               names_to = 'parameter') %>%
  group_by(rowid, parameter) %>%
  summarise(kurtosis = kurtosis(value) - 3,
            range = diff(range(value))) %>%
  ungroup() %>%
  select(rowid, parameter, kurtosis, range)

overall_struq_kurtosis %>% filter_all(any_vars(is.nan(.))) # 138 Ab/DP combinations where DP values are constant--> infinity


overall_struc_sensitivity = mutated_struc_dp_dist %>%
  filter(mutation_region != 'none') %>%
  pivot_longer(-c(rowid, mutation_type, mutation_region, isotype, species),
               values_to = "distance",
               names_to = "parameter") %>%
  group_by(rowid, parameter) %>%
  summarise(ASD = median(distance)) %>%
  ungroup() %>% 
  mutate(parameter = case_when(parameter == "AbStruc_mean_interaction_distance_dist" ~ 'AbStruc_mean_interaction_distance',
                               T ~ gsub("_dist","",parameter))) %>%
  left_join(overall_struq_kurtosis)

# Sensitivity by region
region_seq_kurtosis = mutated_seq_dp %>%
  filter(mutation_region != 'none') %>%
  scale_data() %>%
  pivot_longer(-c(-starts_with('Ab')),
               names_to = 'parameter') %>%
  group_by(rowid, parameter, mutation_region) %>%
  summarise(kurtosis = kurtosis(value) - 3,
            range = diff(range(value))) %>%
  ungroup() %>%
  select(rowid, mutation_region, parameter, kurtosis, range)

region_seq_kurtosis %>% filter_all(any_vars(is.nan(.))) # 1930 infinity due to Ab/DP/region combinations with constant values

region_seq_sensitivity = mutated_seq_dp_dist %>%
  filter(mutation_region != 'none') %>%
  pivot_longer(-c(rowid, mutation_type, mutation_region, isotype, species),
               values_to = "distance",
               names_to = "parameter") %>%
  group_by(rowid, parameter, mutation_region) %>%
  summarise(ASD = mean(distance)) %>%
  ungroup() %>%
  mutate(parameter = case_when(parameter == "AbStruc_mean_interaction_distance_dist" ~ 'AbStruc_mean_interaction_distance',
                               T ~ gsub("_dist","",parameter))) %>%
  left_join(region_seq_kurtosis)

region_struc_kurtosis = mutated_struc_dp %>%
  filter(mutation_region != 'none') %>%
  scale_data() %>%
  pivot_longer(-c(-starts_with('Ab')),
               names_to = 'parameter') %>%
  group_by(rowid, parameter, mutation_region) %>%
  summarise(kurtosis = kurtosis(value) - 3,
            range = diff(range(value))) %>%
  ungroup() %>%
  select(rowid, mutation_region, parameter, kurtosis, range)


region_struc_kurtosis %>% filter_all(any_vars(is.nan(.))) # 492 infinity due to Ab/DP/region combinations with constant values

region_struc_sensitivity = mutated_struc_dp_dist %>%
  filter(mutation_region != 'none') %>%
  pivot_longer(-c(rowid, mutation_type, mutation_region, isotype, species),
               values_to = "distance",
               names_to = "parameter") %>%
  group_by(rowid, parameter, mutation_region) %>%
  summarise(ASD = mean(distance)) %>%
  ungroup() %>%
  mutate(parameter = case_when(parameter == "AbStruc_mean_interaction_distance_dist" ~ 'AbStruc_mean_interaction_distance',
                               T ~ gsub("_dist","",parameter))) %>%
  left_join(region_struc_kurtosis)


# Sensitivity by isotype
isotype_seq_kurtosis = mutated_seq_dp %>%
  filter(mutation_region != 'none') %>%
  scale_data() %>%
  pivot_longer(-c(-starts_with('Ab')),
               names_to = 'parameter') %>%
  group_by(rowid, parameter, isotype) %>%
  summarise(kurtosis = kurtosis(value) - 3,
            range = diff(range(value))) %>%
  ungroup() %>%
  select(rowid, isotype, parameter, kurtosis, range)

isotype_seq_sensitivity = mutated_seq_dp_dist %>%
  filter(mutation_region != 'none') %>%
  pivot_longer(-c(rowid, mutation_type, mutation_region, isotype, species),
               values_to = "distance",
               names_to = "parameter") %>%
  group_by(rowid, parameter, isotype) %>%
  summarise(ASD = mean(distance)) %>%
  ungroup() %>%
  mutate(parameter = case_when(parameter == "AbStruc_mean_interaction_distance_dist" ~ 'AbStruc_mean_interaction_distance',
                               T ~ gsub("_dist","",parameter))) %>%
  left_join(isotype_seq_kurtosis)


isotype_struc_kurtosis = mutated_struc_dp %>%
  filter(mutation_region != 'none') %>%
  scale_data() %>%
  pivot_longer(-c(-starts_with('Ab')),
               names_to = 'parameter') %>%
  group_by(rowid, parameter, isotype) %>%
  summarise(kurtosis = kurtosis(value) - 3,
            range = diff(range(value))) %>%
  ungroup() %>% 
  select(rowid, isotype, parameter, kurtosis, range)

isotype_struc_kurtosis %>% filter_all(any_vars(is.na(.))) # 138 infinities due to Ab/DP/isotype combinations with constant values




isotype_struc_sensitivity = mutated_struc_dp_dist %>%
  filter(mutation_region != 'none') %>%
  pivot_longer(-c(rowid, mutation_type, mutation_region, isotype, species),
               values_to = "distance",
               names_to = "parameter") %>%
  group_by(rowid, parameter, isotype) %>%
  summarise(ASD = mean(distance)) %>%
  ungroup() %>%
  mutate(parameter = case_when(parameter == "AbStruc_mean_interaction_distance_dist" ~ 'AbStruc_mean_interaction_distance',
                               T ~ gsub("_dist","",parameter))) %>% 
  left_join(isotype_struc_kurtosis)


# Sensitivity by species

species_seq_kurtosis = mutated_seq_dp %>%
  filter(mutation_region != 'none') %>%
  scale_data() %>%
  pivot_longer(-c(-starts_with('Ab')),
               names_to = 'parameter') %>%
  group_by(rowid, parameter, species) %>%
  summarise(kurtosis = kurtosis(value) - 3,
            range = diff(range(value))) %>%
  ungroup() %>%
  select(rowid, species, parameter, kurtosis, range)

species_seq_sensitivity = mutated_seq_dp_dist %>%
  filter(mutation_region != 'none') %>%
  pivot_longer(-c(rowid, mutation_type, mutation_region, isotype, species),
               values_to = "distance",
               names_to = "parameter") %>%
  group_by(rowid, parameter, species) %>%
  summarise(ASD = mean(distance)) %>%
  ungroup() %>%
  mutate(parameter = case_when(parameter == "AbStruc_mean_interaction_distance_dist" ~ 'AbStruc_mean_interaction_distance',
                               T ~ gsub("_dist","",parameter))) %>%
  left_join(species_seq_kurtosis)


species_struc_kurtosis = mutated_struc_dp %>%
  filter(mutation_region != 'none') %>%
  scale_data() %>%
  pivot_longer(-c(-starts_with('Ab')),
               names_to = 'parameter') %>%
  group_by(rowid, parameter, species) %>%
  summarise(kurtosis = kurtosis(value) - 3,
            range = diff(range(value))) %>%
  ungroup() %>%
  select(rowid, species, parameter, kurtosis, range)

species_struc_kurtosis %>% filter_all(any_vars(is.na(.))) # 138 infinities due to Ab/DP/species combination with constant values

species_struc_sensitivity = mutated_struc_dp_dist %>%
  filter(mutation_region != 'none') %>%
  pivot_longer(-c(rowid, mutation_type, mutation_region, isotype, species),
               values_to = "distance",
               names_to = "parameter") %>%
  group_by(rowid, parameter, species) %>%
  summarise(ASD = mean(distance)) %>%
  ungroup() %>%
  mutate(parameter = case_when(parameter == "AbStruc_mean_interaction_distance_dist" ~ 'AbStruc_mean_interaction_distance',
                               T ~ gsub("_dist","",parameter))) %>%
  left_join(species_struc_kurtosis)


save(overall_seq_sensitivity,
     overall_struc_sensitivity,
     region_seq_sensitivity,
     region_struc_sensitivity,
     isotype_seq_sensitivity,
     isotype_struc_sensitivity,
     species_seq_sensitivity,
     species_struc_sensitivity,
     file = 'data/human_heavy_sensmetrics.RData')
