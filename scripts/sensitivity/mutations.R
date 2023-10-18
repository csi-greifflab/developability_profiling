library(tidyverse)
library(readr)

df = read_csv("data/developability_data_basic_genes.csv") %>% 
  select(-...1)

rowids_sample = readRDS('data/sensitivity_sample_full_rowid.RDS')

get_mutations = function(df, row_id){
  row = df %>% filter(rowid == row_id)
  number_of_mutations = 19*nchar(row$aaSeqAbChain)
  row_number = 1
  
  #get region range
  fr1_range = str_locate(row$aaSeqAbChain, row$aaSeqFR1)
  fr2_range = str_locate(row$aaSeqAbChain, row$aaSeqFR2)
  fr3_range = str_locate(row$aaSeqAbChain, row$aaSeqFR3)
  fr4_range = str_locate(row$aaSeqAbChain, row$aaSeqFR4)
  cdr1_range = str_locate(row$aaSeqAbChain, row$aaSeqCDR1)
  cdr2_range = str_locate(row$aaSeqAbChain, row$aaSeqCDR2)
  cdr3_range = str_locate(row$aaSeqAbChain, row$aaSeqCDR3)
  
  #create empty new data frame
  new_df = data.frame(matrix(ncol = length(colnames(row)), nrow = number_of_mutations + 1))
  colnames(new_df) = colnames(df)
  new_df = new_df %>% mutate(mutation_type = "")
  new_df[1,] = row[1,]
  new_df[1,'mutation_type'] = "none"
  
  #list of amino acids
  amino_acids = c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
  
  #iterates over full sequence 
  for(pos in 1: nchar(row$aaSeqAbChain)){
    #iterates over all amino acids
    for(aa in amino_acids){
      if(aa != substr(row$aaSeqAbChain, pos, pos)){
        #Set new index
        row_number = row_number + 1
        
        #add mutation type and mutated full sequence
        new_df$mutation_type[row_number] = paste(substr(row$aaSeqAbChain, pos, pos),pos,aa,sep="")
        mutated_seq = row$aaSeqAbChain
        substr(mutated_seq, pos, pos) = aa
        new_df$aaSeqAbChain[row_number] = mutated_seq
      }
    }
  }
  
  #clean up new_df
  new_df = new_df %>% 
    #copy non sequence values do all other rows
    mutate(isotype = pluck(new_df, "isotype",1)) %>%
    mutate(database = pluck(new_df, "database",1)) %>%
    mutate(species = pluck(new_df, "species",1)) %>%
    mutate(chain = pluck(new_df, "chain",1)) %>%
    mutate(v_gene = pluck(new_df, "v_gene",1)) %>%
    mutate(j_gene = pluck(new_df, "j_gene",1)) %>%
    mutate(rowid = pluck(new_df, "rowid",1)) %>%
    #divide full sequence into subregions
    mutate(aaSeqFR1 = substr(new_df$aaSeqAbChain, fr1_range[1], fr1_range[2])) %>%
    mutate(aaSeqFR2 = substr(new_df$aaSeqAbChain, fr2_range[1], fr2_range[2])) %>%
    mutate(aaSeqFR3 = substr(new_df$aaSeqAbChain, fr3_range[1], fr3_range[2])) %>%
    mutate(aaSeqFR4 = substr(new_df$aaSeqAbChain, fr4_range[1], fr4_range[2])) %>%
    mutate(aaSeqCDR1 = substr(new_df$aaSeqAbChain, cdr1_range[1], cdr1_range[2])) %>%
    mutate(aaSeqCDR2 = substr(new_df$aaSeqAbChain, cdr2_range[1], cdr2_range[2])) %>%
    mutate(aaSeqCDR3 = substr(new_df$aaSeqAbChain, cdr3_range[1], cdr3_range[2]))
  return(new_df)
}

mutated_datalist = list()

for (id in rowids_sample) {
  mutated_datalist[[as.character(id)]] = get_mutations(df, id)
}

mutated_dataset_full = rbindlist(mutated_datalist) %>% as_tibble()
saveRDS(mutated_dataset_full, '/storage/jahnz/R/developability/data/sensitivity/mutated_dataset_full.RDS')


rowids_subsample = readRDS('/storage/jahnz/R/developability/data/sensitivity/sensitivity_sample_subsample_rowid.RDS')

mutated_dataset_sampled = mutated_dataset_full %>% filter(rowid %in% rowids_subsample)

write.csv(mutated_dataset_sampled, file = '/storage/jahnz/R/developability/data/sensitivity/mutated_dataset_sampled.csv')






#Calculate Parameters

library(parallel)
library(stringr)
library(Peptides)

calculated_parameters = mutated_dataset_full
cl <- makeCluster(16)

#####

#get sequence lengths and molecular weights:

calculated_parameters = calculated_parameters %>%
  mutate(cdr1_length = nchar(aaSeqCDR1),
         cdr2_length = nchar(aaSeqCDR2),
         cdr3_length = nchar(aaSeqCDR3),
         AbChain_length = nchar(aaSeqAbChain),
         fr1_length = nchar(aaSeqFR1),
         fr2_length = nchar(aaSeqFR2),
         fr3_length = nchar(aaSeqFR3),
         fr4_length = nchar(aaSeqFR4))

calculated_parameters = calculated_parameters %>%
  mutate(AbChain_mw = mw(aaSeqAbChain, monoisotopic = FALSE),
         fr1_mw = mw(aaSeqFR1, monoisotopic = FALSE),
         fr2_mw = mw(aaSeqFR2, monoisotopic = FALSE),
         fr3_mw = mw(aaSeqFR3, monoisotopic = FALSE),
         fr4_mw = mw(aaSeqFR4, monoisotopic = FALSE),
         cdr1_mw = mw(aaSeqCDR1, monoisotopic = FALSE),
         cdr2_mw = mw(aaSeqCDR2, monoisotopic = FALSE),
         cdr3_mw = mw(aaSeqCDR3, monoisotopic = FALSE)) 

calculated_parameters = calculated_parameters %>%
  mutate(fr1_av_resdiue_weight = fr1_mw/(fr1_length - 1),
         fr2_av_resdiue_weight = fr2_mw/(fr2_length - 1),
         fr3_av_resdiue_weight = fr3_mw/(fr3_length - 1),
         fr4_av_resdiue_weight = fr4_mw/(fr4_length - 1),
         cdr1_av_resdiue_weight = cdr1_mw/(cdr1_length - 1),
         cdr2_av_resdiue_weight = cdr2_mw/(cdr2_length - 1),
         cdr3_av_resdiue_weight = cdr3_mw/(cdr3_length - 1),
         AbChain_av_resdiue_weight = AbChain_mw/(AbChain_length - 1))
#####
#Instaindex

#add the hydrophobicity to the basic database
print("start Instaindex calculation")

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

calculated_parameters = calculated_parameters %>%
  mutate(AbChain_instaindex = unlist(clusterApply(cl, calculated_parameters$aaSeqAbChain, fun=get_instaindex)),
         cdr1_instaindex = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR1, fun=get_instaindex)),
         cdr2_instaindex = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR2, fun=get_instaindex)),
         cdr3_instaindex = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR3, fun=get_instaindex)),
         fr1_instaindex = unlist(clusterApply(cl, calculated_parameters$aaSeqFR1, fun=get_instaindex)),
         fr2_instaindex = unlist(clusterApply(cl, calculated_parameters$aaSeqFR2, fun=get_instaindex)),
         fr3_instaindex = unlist(clusterApply(cl, calculated_parameters$aaSeqFR3, fun=get_instaindex)),
         fr4_instaindex = unlist(clusterApply(cl, calculated_parameters$aaSeqFR4, fun=get_instaindex)))

get_aliphindex = function(seq){
  library(parallel)
  library(stringr)
  library(Peptides)
  aIndex(seq)
}

calculated_parameters = calculated_parameters %>%
  mutate(AbChain_aliphindex = unlist(clusterApply(cl, calculated_parameters$aaSeqAbChain, fun=get_aliphindex)),
         cdr1_aliphindex = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR1, fun=get_aliphindex)),
         cdr2_aliphindex = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR2, fun=get_aliphindex)),
         cdr3_aliphindex = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR3, fun=get_aliphindex)),
         fr1_aliphindex = unlist(clusterApply(cl, calculated_parameters$aaSeqFR1, fun=get_aliphindex)),
         fr2_aliphindex = unlist(clusterApply(cl, calculated_parameters$aaSeqFR2, fun=get_aliphindex)),
         fr3_aliphindex = unlist(clusterApply(cl, calculated_parameters$aaSeqFR3, fun=get_aliphindex)),
         fr4_aliphindex = unlist(clusterApply(cl, calculated_parameters$aaSeqFR4, fun=get_aliphindex)))

print("Instaindex done")

#####
#Hydrophobicity

#add the hydrophobicity to the basic database
print("start Hydrophobicity calculations")

get_hydrophobicity = function(seq){
  library(parallel)
  library(stringr)
  library(Peptides)
  hydrophobicity(seq, scale = "Eisenberg")
}

calculated_parameters = calculated_parameters %>%
  mutate(AbChain_hydrophobicity = unlist(clusterApply(cl, calculated_parameters$aaSeqAbChain, fun=get_hydrophobicity)),
         cdr1_hydrophobicity = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR1, fun=get_hydrophobicity)),
         cdr2_hydrophobicity = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR2, fun=get_hydrophobicity)),
         cdr3_hydrophobicity = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR3, fun=get_hydrophobicity)),
         fr1_hydrophobicity = unlist(clusterApply(cl, calculated_parameters$aaSeqFR1, fun=get_hydrophobicity)),
         fr2_hydrophobicity = unlist(clusterApply(cl, calculated_parameters$aaSeqFR2, fun=get_hydrophobicity)),
         fr3_hydrophobicity = unlist(clusterApply(cl, calculated_parameters$aaSeqFR3, fun=get_hydrophobicity)),
         fr4_hydrophobicity = unlist(clusterApply(cl, calculated_parameters$aaSeqFR4, fun=get_hydrophobicity)))

print("Hydrophobicity done")
#####
#Hydrophobic moment
print("Start hydrophobic moment calculation")

get_hmom = function(seq_input){
  library(parallel)
  library(stringr)
  library(Peptides)
  hmoment(seq = seq_input, angle = 160, window = 10)
}

calculated_parameters = calculated_parameters %>%
  mutate(AbChain_hmom = unlist(clusterApply(cl, calculated_parameters$aaSeqAbChain, fun=get_hmom)))

print("Hydrophobic moment done")

#####
#Charge

charge_intervals = as.character(c(1:14))


get_charge = function(seq_input, pH_input){
  library(parallel)
  library(stringr)
  library(Peptides)
  charge_intervals = as.character(c(1:14))
  charge(seq_input, pH = pH_input , pKscale = "Lehninger")
}

for (i in 1:length(charge_intervals)){
  charge_intervals = as.character(c(1:14))
  print(charge_intervals[i])
  pH_input = as.numeric(charge_intervals[i])
  calculated_parameters = calculated_parameters %>%
    #sample_n(2) %>%
    mutate(!!paste0("AbChain_",charge_intervals[i],"_charge") := unlist(clusterMap(cl, get_charge, seq_input = calculated_parameters$aaSeqAbChain, pH_input = pH_input)),
           !!paste0("cdr1_",charge_intervals[i],"_charge") := unlist(clusterMap(cl, get_charge, seq_input = calculated_parameters$aaSeqCDR1, pH_input = pH_input)),
           !!paste0("cdr2_",charge_intervals[i],"_charge") := unlist(clusterMap(cl, get_charge, seq_input = calculated_parameters$aaSeqCDR2, pH_input = pH_input)),
           !!paste0("cdr3_",charge_intervals[i],"_charge") := unlist(clusterMap(cl, get_charge, seq_input = calculated_parameters$aaSeqCDR3, pH_input = pH_input)),
           !!paste0("fr1_",charge_intervals[i],"_charge") := unlist(clusterMap(cl, get_charge, seq_input = calculated_parameters$aaSeqFR1, pH_input = pH_input)),
           !!paste0("fr2_",charge_intervals[i],"_charge") := unlist(clusterMap(cl, get_charge, seq_input = calculated_parameters$aaSeqFR2, pH_input = pH_input)),
           !!paste0("fr3_",charge_intervals[i],"_charge") := unlist(clusterMap(cl, get_charge, seq_input = calculated_parameters$aaSeqFR3, pH_input = pH_input)),
           !!paste0("fr4_",charge_intervals[i],"_charge") := unlist(clusterMap(cl, get_charge, seq_input = calculated_parameters$aaSeqFR4, pH_input = pH_input)))
}

#####
#Aromaticity

print("Start aromaticity calculation")
get_Aromaticity = function(input){
  library(parallel)
  library(stringr)
  library(Peptides)
  output = ((str_count(input, "F") + str_count(input,"H") + str_count(input, "W") + str_count(input,"Y"))/nchar(input))*100
}


calculated_parameters = calculated_parameters %>% 
  mutate(fr1_aromaticity = unlist(clusterApply(cl, calculated_parameters$aaSeqFR1, fun=get_Aromaticity)),
         fr2_aromaticity = unlist(clusterApply(cl, calculated_parameters$aaSeqFR2, fun=get_Aromaticity)),
         fr3_aromaticity = unlist(clusterApply(cl, calculated_parameters$aaSeqFR3, fun=get_Aromaticity)),
         fr4_aromaticity = unlist(clusterApply(cl, calculated_parameters$aaSeqFR4, fun=get_Aromaticity))) 

calculated_parameters = calculated_parameters %>% 
  mutate(cdr1_aromaticity = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR1, fun=get_Aromaticity)),
         cdr2_aromaticity = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR2, fun=get_Aromaticity)),
         cdr3_aromaticity = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR3, fun=get_Aromaticity)))
print("Aromaticity calculation done")

#####
#Absorbance molecular extinction coefficient

print("Start molecular extinction coefficient calculation")
#This script is togenerate light absorbance info (ext coefficients)

get_mol_extcoef = function(input){
  library(parallel)
  library(stringr)
  library(Peptides)
  if(str_count(input,"C") %% 2 == 0){
    (str_count(input, "Y")*1490) + (str_count(input,"W")*5500) + ((str_count(input,"C")/2)*125)
  }
  else{
    (str_count(input, "Y")*1490) + (str_count(input,"W")*5500) + (((str_count(input,"C") -1)  /2)*125)
  }
} 


calculated_parameters = calculated_parameters %>%
  mutate(AbChain_molextcoef = unlist(clusterApply(cl, calculated_parameters$aaSeqAbChain, fun=get_mol_extcoef)),
         cdr1_molextcoef = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR1, fun=get_mol_extcoef)),
         cdr2_molextcoef = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR2, fun=get_mol_extcoef)),
         cdr3_molextcoef = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR3, fun=get_mol_extcoef)),
         fr1_molextcoef = unlist(clusterApply(cl, calculated_parameters$aaSeqFR1, fun=get_mol_extcoef)),
         fr2_molextcoef = unlist(clusterApply(cl, calculated_parameters$aaSeqFR2, fun=get_mol_extcoef)),
         fr3_molextcoef = unlist(clusterApply(cl, calculated_parameters$aaSeqFR3, fun=get_mol_extcoef)),
         fr4_molextcoef = unlist(clusterApply(cl, calculated_parameters$aaSeqFR4, fun=get_mol_extcoef)))

get_mol_extcoef_cystine_bridges = function(input){
  library(parallel)
  library(stringr)
  library(Peptides)
  if(str_count(input,"C") %% 2 == 0){
    ((str_count(input,"C")/2)*125)
  }
  else{
    (((str_count(input,"C") -1)  /2)*125)
  }
} 


calculated_parameters = calculated_parameters %>%
  mutate(AbChain_cysbridges_molextcoef = unlist(clusterApply(cl, calculated_parameters$aaSeqAbChain, fun=get_mol_extcoef_cystine_bridges)),
         cdr1_cysbridges_molextcoef = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR1, fun=get_mol_extcoef_cystine_bridges)),
         cdr2_cysbridges_molextcoef = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR2, fun=get_mol_extcoef_cystine_bridges)),
         cdr3_cysbridges_molextcoef = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR3, fun=get_mol_extcoef_cystine_bridges)),
         fr1_cysbridges_molextcoef = unlist(clusterApply(cl, calculated_parameters$aaSeqFR1, fun=get_mol_extcoef_cystine_bridges)),
         fr2_cysbridges_molextcoef = unlist(clusterApply(cl, calculated_parameters$aaSeqFR2, fun=get_mol_extcoef_cystine_bridges)),
         fr3_cysbridges_molextcoef = unlist(clusterApply(cl, calculated_parameters$aaSeqFR3, fun=get_mol_extcoef_cystine_bridges)),
         fr4_cysbridges_molextcoef = unlist(clusterApply(cl, calculated_parameters$aaSeqFR4, fun=get_mol_extcoef_cystine_bridges)))

print("Molecular extinction coefficient done")

#####
#Absorbance Percentage

print("Start absorbance percentage calculation")

calculated_parameters = mutated_dataset_full %>%
  select(c(rowid, "aaSeqAbChain", "aaSeqFR1", "aaSeqCDR1", "aaSeqFR2", "aaSeqCDR2", "aaSeqFR3", "aaSeqCDR3", "aaSeqFR4", "isotype", "database", "species", "chain", ends_with("mw"))) %>%
  left_join(calculated_parameters) %>%
  mutate(AbChain_percentextcoef = (AbChain_molextcoef*10)/AbChain_mw,
         cdr1_percentextcoef = (cdr1_molextcoef*10)/cdr1_mw,
         cdr2_percentextcoef = (cdr2_molextcoef*10)/cdr2_mw,
         cdr3_percentextcoef = (cdr3_molextcoef*10)/cdr3_mw,
         fr1_percentextcoef = (fr1_molextcoef*10)/fr1_mw,
         fr2_percentextcoef = (fr2_molextcoef*10)/fr2_mw,
         fr3_percentextcoef = (fr3_molextcoef*10)/fr3_mw,
         fr4_percentextcoef = (fr4_molextcoef*10)/fr4_mw,
         AbChain_cysbridges_percentextcoef = (AbChain_cysbridges_molextcoef*10)/AbChain_mw,
         cdr1_cysbridges_percentextcoef = (cdr1_cysbridges_molextcoef*10)/cdr1_mw,
         cdr2_cysbridges_percentextcoef = (cdr2_cysbridges_molextcoef*10)/cdr2_mw,
         cdr3_cysbridges_percentextcoef = (cdr3_cysbridges_molextcoef*10)/cdr3_mw,
         fr1_cysbridges_percentextcoef = (fr1_cysbridges_molextcoef*10)/fr1_mw,
         fr2_cysbridges_percentextcoef = (fr2_cysbridges_molextcoef*10)/fr2_mw,
         fr3_cysbridges_percentextcoef = (fr3_cysbridges_molextcoef*10)/fr3_mw,
         fr4_cysbridges_percentextcoef = (fr4_cysbridges_molextcoef*10)/fr4_mw)

print("Absorbance percentage calculation done")
#####
#Isoelectric point calculation

print("Absorbance percentage calculation done")
get_pI = function(seq_input){
  library(Peptides)
  library(tidyverse)
  library(parallel)
  pI(seq_input, pKscale = "Lehninger")
}  

calculated_parameters = calculated_parameters %>%
  mutate(AbChain_pI = unlist(clusterApply(cl, calculated_parameters$aaSeqAbChain, fun=get_pI)),
         cdr1_pI = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR1, fun=get_pI)),
         cdr2_pI = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR2, fun=get_pI)),
         cdr3_pI = unlist(clusterApply(cl, calculated_parameters$aaSeqCDR3, fun=get_pI)),
         fr1_pI = unlist(clusterApply(cl, calculated_parameters$aaSeqFR1, fun=get_pI)),
         fr2_pI = unlist(clusterApply(cl, calculated_parameters$aaSeqFR2, fun=get_pI)),
         fr3_pI = unlist(clusterApply(cl, calculated_parameters$aaSeqFR3, fun=get_pI)),
         fr4_pI = unlist(clusterApply(cl, calculated_parameters$aaSeqFR4, fun=get_pI)))

print("Absorbance percentage calculation done")
print("All calculations finished")

saveRDS(calculated_parameters, "data/parameters_mutations.rds")
