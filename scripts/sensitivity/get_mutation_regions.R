library(tidyverse)
library(data.table)
library(stringr)


AbChain_mutated_full_developability = read_csv("data/AbChain_mutated_full_developability.csv")

mutation_regions = AbChain_mutated_full_developability %>% select(c('rowid', 'aaSeqFR1', 'aaSeqCDR1', 'aaSeqFR2', 'aaSeqCDR2', 'aaSeqFR3', 'aaSeqCDR3', 'aaSeqFR4', 'mutation_type'))

# Isolate mutation position 
mutation_regions = mutation_regions %>% mutate(mutation_position = case_when(mutation_type == 'none' ~ 'none',
                                                                             TRUE ~ str_sub(mutation_type,2,-2)))
                                               


# Convert region sequences to region length and drop sequences
mutation_regions_reduced = mutation_regions %>% 
  mutate(FR1_length = nchar(aaSeqFR1)) %>%
  mutate(CDR1_length = nchar(aaSeqCDR1)) %>%
  mutate(FR2_length = nchar(aaSeqFR2)) %>%
  mutate(CDR2_length = nchar(aaSeqCDR2)) %>%
  mutate(FR3_length = nchar(aaSeqFR3)) %>%
  mutate(CDR3_length = nchar(aaSeqCDR3)) %>%
  mutate(FR4_length = nchar(aaSeqFR4)) %>%
  select(c(rowid, mutation_type, mutation_position, FR1_length, CDR1_length, FR2_length, CDR2_length, FR3_length, CDR3_length, FR4_length))

# Convert position column to integer  
mutation_regions_reduced = mutation_regions_reduced %>% mutate(mutation_position = as.integer(mutation_position))

# Checks if mutation position is smaller or equal to region length, if not, it compares the next region length + all previous regions and so on.
mutation_regions_reduced = mutation_regions_reduced %>% mutate(mutation_region = case_when(mutation_type == 'none' ~ 'none',
                                                                                           mutation_position <= FR1_length ~ 'FR1',
                                                                                           mutation_position <= FR1_length + CDR1_length ~ 'CDR1',
                                                                                           mutation_position <= FR1_length + CDR1_length + FR2_length ~ 'FR2',
                                                                                           mutation_position <= FR1_length + CDR1_length + FR2_length + CDR2_length ~ 'CDR2',
                                                                                           mutation_position <= FR1_length + CDR1_length + FR2_length + CDR2_length + FR3_length ~ 'FR3',
                                                                                           mutation_position <= FR1_length + CDR1_length + FR2_length + CDR2_length + FR3_length + CDR3_length ~ 'CDR3',
                                                                                           TRUE ~ 'FR4'))

mutation_regions_for_merge = mutation_regions_reduced %>% select(rowid, mutation_type, mutation_region)
colnames(AbChain_mutated_full_developability)
colnames(mutation_regions_for_merge)

AbChain_mutated_full_developability = merge(AbChain_mutated_full_developability,mutation_regions_for_merge, by = c('rowid', 'mutation_type')) %>% as_tibble()
saveRDS(AbChain_mutated_full_developability, "data/AbChain_mutated_full_developability.rds")


