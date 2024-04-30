#Author: Habib Bashour

library(Peptides)
library(tidyverse)
library(parallel)
library(doParallel)
library(foreach)
library(stringr)

#first load/read a dataframe that contains the sequences you are interested in

developability_data = read_csv("~/data/dataframe.csv") #this data should include aaSeqAbChain column that contains the sequence of the variable region (fv) of the antibody.

# Specify number of cores for parallelization
num_cores = 50
cl <- makeCluster(num_cores)

#charge----

charge_intervals = as.character(c(1:14))


for (i in 1:length(charge_intervals)){
  charge_intervals = as.character(c(1:14))
  print(charge_intervals[i])
  pH_input = as.numeric(charge_intervals[i])
  get_charge = function(seq_input, pH_input){
    charge_intervals = as.character(c(1:14))
    library(Peptides)
    library(tidyverse)
    library(parallel)
    charge(seq_input, pH = pH_input , pKscale = "Lehninger")
  }
  developability_data = developability_data %>%
    mutate(!!paste0("AbChain_",charge_intervals[i],"_charge") := unlist(clusterMap(cl, get_charge, seq_input = developability_data$aaSeqAbChain, pH_input = pH_input)))
}

#pI----

get_pI = function(seq_input){
  library(Peptides)
  library(tidyverse)
  library(parallel)
  pI(seq_input, pKscale = "Lehninger")
}  


developability_data = developability_data %>%
  mutate(AbChain_pI = unlist(clusterApply(cl, developability_data$aaSeqAbChain, fun=get_pI)))

#photochemical----

get_mol_extcoef = function(input){
  library(stringr)
  if(str_count(input,"C") %% 2 == 0){
    (str_count(input, "Y")*1490) + (str_count(input,"W")*5500) + ((str_count(input,"C")/2)*125)
  }
  else{
    (str_count(input, "Y")*1490) + (str_count(input,"W")*5500) + (((str_count(input,"C") -1)  /2)*125)
  }
} 

developability_data = developability_data %>%
  mutate(AbChain_molextcoef = unlist(clusterApply(cl, developability_data$aaSeqAbChain, fun=get_mol_extcoef)))


get_mol_extcoef_cystine_bridges = function(input){
  library(stringr)
  if(str_count(input,"C") %% 2 == 0){
    ((str_count(input,"C")/2)*125)
  }
  else{
    (((str_count(input,"C") -1)  /2)*125)
  }
} 

developability_data = developability_data %>%
  mutate(AbChain_cysbridges_molextcoef = unlist(clusterApply(cl, developability_data$aaSeqAbChain, fun=get_mol_extcoef_cystine_bridges)))

developability_data = developability_data %>%
  mutate(AbChain_percentextcoef = (AbChain_molextcoef*10)/AbChain_mw) %>%
  mutate(AbChain_cysbridges_percenextcoef = (AbChain_cysbridges_molextcoef*10)/AbChain_mw)


#indexes---- 
get_instaindex = function(seq){
  library(Peptides)
  library(tidyverse)
  library(parallel)
  if(nchar(seq)<3){
    return(as.double("NA"))
  }
  else {
    return(instaIndex(seq))
  }
}
get_aliphindex = function(seq){
  library(Peptides)
  library(tidyverse)
  library(parallel)
  aIndex(seq)
}


developability_data = developability_data %>%
  mutate(AbChain_instaindex = unlist(clusterApply(cl, developability_data$aaSeqAbChain, fun=get_instaindex)),
         AbChain_aliphindex = unlist(clusterApply(cl, developability_data$aaSeqAbChain, fun=get_aliphindex)))


#hydrophobicity and hmom----

get_hydrophobicity = function(seq){
  library(Peptides)
  library(tidyverse)
  library(parallel)
  hydrophobicity(seq, scale = "Eisenberg")
}

developability_data = developability_data %>%
  mutate(AbChain_hydrophobicity = unlist(clusterApply(cl, developability_data$aaSeqAbChain, fun=get_hydrophobicity)))


get_hmom = function(seq_input){
  library(Peptides)
  library(tidyverse)
  library(parallel)
  hmoment(seq = seq_input, angle = 160, window = 10)
}


developability_data = developability_data %>%
  mutate(AbChain_hmom = unlist(clusterApply(cl, developability_data$aaSeqAbChain, fun=get_hmom)))


#amino acid categorical content----- 

get_Aromaticity = function(input){
  library(tidyverse)
  library(parallel)
  output = ((str_count(input, "F") + str_count(input,"H") + str_count(input, "W") + str_count(input,"Y"))/nchar(input))*100
}

# Tiny (A+C+G+S+T)
get_Tiny = function(input){
  library(parallel)
  library(stringr)
  output = ((str_count(input, "A") + str_count(input,"C") + str_count(input, "G") + str_count(input,"S") + str_count(input,"T") )/nchar(input))*100
}
# Small (A+C+D+G+N+P+S+T+V)

get_Small = function(input){
  library(stringr)
  library(parallel)
  output = ((str_count(input, "A") + str_count(input,"C") + str_count(input, "D") +  str_count(input,"G") + str_count(input,"N") + str_count(input,"P") + str_count(input,"S") + str_count(input,"T") + str_count(input,"V") )/nchar(input))*100
}
# Aliphatic (A+I+L+V)

get_Aliphatic = function(input){
  library(stringr)
  library(parallel)
  output = ((str_count(input, "A") + str_count(input,"I") + str_count(input, "L") + str_count(input,"V") )/nchar(input))*100
}

# Nonpolar (A+C+F+G+I+L+M+P+V+W+Y) 
get_Nonpolar = function(input){
  library(stringr)
  library(parallel)
  output = ((str_count(input, "A") + str_count(input,"C") + str_count(input, "F") + str_count(input,"G") + str_count(input,"I") + 
               str_count(input,"L") + str_count(input,"M") + str_count(input,"P") + str_count(input,"V") + str_count(input,"W") + str_count(input,"Y"))/nchar(input))*100
}

# Polar  (D+E+H+K+N+Q+R+S+T)
get_Polar = function(input){
  library(stringr)
  library(parallel)
  output = ((str_count(input, "D") + str_count(input,"E") + str_count(input, "H") + str_count(input,"K") + str_count(input,"N") + 
               str_count(input,"Q") + str_count(input,"R") + str_count(input,"S") + str_count(input,"T"))/nchar(input))*100
}

# Basic (H+K+R)
get_Basic = function(input){
  library(stringr)
  library(parallel)
  output = ((str_count(input,"H") + str_count(input, "K") + str_count(input,"R"))/nchar(input))*100
}

# Acidic (D+E)
get_Acidic = function(input){
  library(stringr)
  library(parallel)
  output = ((str_count(input,"D") + str_count(input, "E"))/nchar(input))*100
}


developability_data = developability_data %>% 
  mutate(AbChain_aromatic_content = unlist(clusterApply(cl, developability_data$aaSeqAbChain, fun=get_Aromaticity)),
         AbChain_tiny_content = unlist(clusterApply(cl, developability_data$aaSeqAbChain, fun=get_Tiny)),
         AbChain_small_content = unlist(clusterApply(cl, developability_data$aaSeqAbChain, fun=get_Small)),
         AbChain_aliphatic_content = unlist(clusterApply(cl, developability_data$aaSeqAbChain, fun=get_Aliphatic)),
         AbChain_nonpolar_content = unlist(clusterApply(cl, developability_data$aaSeqAbChain, fun=get_Nonpolar)),
         AbChain_polar_content = unlist(clusterApply(cl, developability_data$aaSeqAbChain, fun=get_Polar)),
         AbChain_basic_content = unlist(clusterApply(cl, developability_data$aaSeqAbChain, fun=get_Basic)),
         AbChain_acidic_content = unlist(clusterApply(cl, developability_data$aaSeqAbChain, fun=get_Acidic)))


#immunogenicity----- 
library(janitor)

#immunogenicity was calculated on 10% only of the native dataset, 105 of PAD. 10% of kymouse, all mAbs (from TheraSAbDAb) 

# (1) produce the input (FASTA):

fastafile = "~/Ab_developability_project/written_fastas/input_file.fasta"
append=FALSE
for(i in 1:nrow(developability_data)){
  if(i > 1){append=TRUE}
  write(paste0(">",developability_data[i,]$isotype, "|" ,developability_data[i,]$species, "|", developability_data[i,]$row_id), # #sequence
        file=fastafile,append=append)
  write(developability_data[i,]$aaSeqAbChain,
        file=fastafile,append=TRUE)
}

# (2) run the immunogenicity command on Terminal using the script (netMHCIIpan.txt)

# (3) read the output (.xls) using these functions 

get_colnames_alleles <- function(file) {
  first_row <- read_lines(file, n_max = 1) %>%
    str_split("\t") %>%
    unlist()
  
  for(i in 1:length(first_row)) {
    if (first_row[i] != "") {
      if ((i %% 2) != 0) {
        first_row[i + 1] = str_c(first_row[i], ":")
        first_row[i]     = str_c(first_row[i], ":")
        
      }
    }
  }
  
  second_row <- read_lines(file, skip = 1, n_max = 1) %>%
    str_split("\t") %>%
    unlist()
  
  str_c(first_row, second_row)
}
read_mhc2 <- function(file) {
  read_tsv(file, skip = 2,
           col_names = get_colnames_alleles(file)) %>%
    pivot_longer(-c(Pos:Target, Ave, NB), names_to = c("allele", "statistic"), names_sep = ":")
}

immuno_data = read_mhc2(netMHCIIpan_output.xls)

get_immunogenicity_measures = function(data_input){
  min_rank <- data_input %>%
    filter(statistic =="Rank") %>%
    group_by(ID) %>%
    #muate(strength())
    slice_min(value) %>%
    rename(min_rank = value) %>%
    distinct(ID, .keep_all = TRUE)
  
  strength_count <- data_input %>%
    filter(statistic =="Rank") %>%
    mutate(strength = case_when(value < 2 ~ "strong",
                                value < 10 ~ "weak",
                                TRUE ~ "NA")) %>%
    group_by(ID, strength) %>%
    summarise(count = n()) %>%
    pivot_wider(names_from = "strength",values_from = "count") %>%
    select(-`NA`) %>%
    replace_na(list(weak = 0, strong = 0)) %>%
    rename(num_weak_binders = weak, num_strong_binder = strong)
  
  average_total_rank <- data_input %>%
    filter(statistic =="Rank") %>%
    group_by(Peptide, ID) %>%
    summarise(peptide_average_rank = mean(value)) %>%
    group_by(ID) %>%
    summarise(full_average_rank = mean(peptide_average_rank))
  
  result = left_join(min_rank, strength_count) %>%
    left_join(average_total_rank) %>%
    rename(maximal_immuno_allele = allele) %>%
    select(-c(NB, Target, Ave)) %>%
    rename(maximal_immuno_position = Pos)
  return(result)
  
}

immuno_data = get_immunogenicity_measures(immuno_data) %>%
  separate(ID, c("isotype","species","row_id"), sep="\\|") %>%
  mutate(rowid = as.integer(rowid))



developability_data = left_join(developability_data, immuno_data)

# Get antibody region of immunogenic peptides

get_regions = function (df, id) {
  filtered = df %>% filter(rowid == id)
  full_seq = pull(filtered, aaSeqAbChain)
  fr1_seq = pull(filtered, aaSeqFR1)
  fr2_seq = pull(filtered, aaSeqFR2)
  fr3_seq = pull(filtered, aaSeqFR3)
  fr4_seq = pull(filtered, aaSeqFR4)
  cdr1_seq = pull(filtered, aaSeqCDR1)
  cdr2_seq = pull(filtered, aaSeqCDR2)
  cdr3_seq = pull(filtered, aaSeqCDR3)
  
  fr1_pos = str_locate(full_seq, fr1_seq)[,2]
  fr2_pos = str_locate(full_seq, fr2_seq)[,2]
  fr3_pos = str_locate(full_seq, fr3_seq)[,2]
  fr4_pos = str_locate(full_seq, fr4_seq)[,2]
  cdr1_pos = str_locate(full_seq, cdr1_seq)[,2]
  cdr2_pos = str_locate(full_seq, cdr2_seq)[,2]
  cdr3_pos = str_locate(full_seq, cdr3_seq)[,2]
  positions = list(id, fr1_pos, cdr1_pos, fr2_pos, cdr2_pos, fr3_pos, cdr3_pos, fr4_pos)
  return(positions)
}

get_peptide_position = function (region_pos) {
  # Calculate peptide start and end regions and numbers of regions spanned by peptide
  new_df = foreach (row = c(1:length(region_pos$rowid)), .combine = 'rbind', .init = tibble()) %dopar% {
    end_not_found = TRUE
    # Region 1 = FR1, Region 2 = CDR1, Region 3 = FR2, etc....
    region = 1
    while (end_not_found) {
      # Increment region by one until it is matched (end_not_found == FALSE)
      region = region + 1
      # if peptide start region is smaller or equal to current region
      if (region_pos[row,]$peptide_start <= region_pos[row,][region]) {
        
        #if peptide end region is smaller or equal to current region
        if (region_pos[row,]$peptide_end <= region_pos[row,][region]) {
          
          #convert region number to region name and assign
          peptide_start_region = colnames(region_pos)[region]
          peptide_end_region = colnames(region_pos)[region]
          peptide_regions_spanned = 1
          end_not_found = FALSE
          region = 2
        }
        else {
          for (end_region in (region+1):8) {
            if (region_pos[row,]$peptide_end <= region_pos[row,][end_region]) {
              peptide_regions_spanned = end_region - region + 1
              peptide_start_region = colnames(region_pos)[region]
              peptide_end_region = colnames(region_pos)[end_region]
              end_not_found = FALSE
              region = 2
              break
            }
          }
        }
      }
    }
    
    frame = data.frame(as.double(region_pos[row, 1]),
                       peptide_start_region,
                       peptide_end_region,
                       peptide_regions_spanned)
    return(frame)
  }
  
  # Rename output data frame
  colnames(new_df) = c("rowid", 
                       "peptide_start_region", 
                       "peptide_end_region", 
                       "peptide_regions_spanned")
  
  return(new_df)
}

get_reg_and_pept = function (df) {
  temp_list = list(get_regions(df, df$rowid))
  temp_frame = data.frame(temp_list) %>% as_tibble()
  colnames(temp_frame) = c("rowid", "fr1", "cdr1", "fr2", "cdr2", "fr3", "cdr3", "fr4")
  
  peptide = pull(df, Peptide)
  rowid = pull(df, rowid)
  peptide_start = pull(df, AbChain_maximal_immuno_position)
  peptide_end = peptide_start + nchar(peptide) - 1
  peptide_df = tibble(data.frame(rowid, peptide_start, peptide_end))
  
  temp_frame = left_join(temp_frame, peptide_df)
  
  return(get_peptide_position(temp_frame))
}

immunogenicity_regions_span = developability_data %>%
  filter(!is.na(AbChain_maximal_immuno_position)) %>% 
  get_reg_and_pept() %>% 
  rename(rowid = rowid,
         AbChain_immunopeptide_regions_span = peptide_regions_spanned)

developability_data = developability_data %>% left_join(immunogenicity_regions_span)

#solubility----
#(1) produce the input file as above (line 210) 

#(2) run the solubility command on Terminal using the script (SoluProt.txt)

#(3) read the output (.csv)
solubility_data = read_csv("solubility_result.csv") %>%
  select(-runtime_id) %>%
  separate(fa_id, c("isotype","species","row_id"), sep  = "\\|") %>%
  rename(AbChain_solubility = soluble)


developability_data = left_join(developability_data, solubility_data)




