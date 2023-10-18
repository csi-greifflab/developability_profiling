library(tidyverse)
library(parallel)
library(stringr)
library(Peptides)

humanised_mice <- read_csv("data/humanised_mice_sampled210K.csv")
output_file = "data/AbChain_developability_humanized_mice"

calculated_parameters = humanised_mice
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

calculated_parameters = calculated_parameters %>%
  select(c("aaSeqAbChain", "aaSeqFR1", "aaSeqCDR1", "aaSeqFR2", "aaSeqCDR2", "aaSeqFR3", "aaSeqCDR3", "aaSeqFR4", "isotype", "database", "species", "chain", ends_with("mw"))) %>%
  left_join(calculated_parameters) %>%
  mutate(AbChain_percentextcoef = (AbChain_molextcoef*10)/AbChain_mw,
         cdr1_percentextcoef = (cdr1_molextcoef*10)/cdr1_mw,
         cdr2_percentextcoef = (cdr2_molextcoef*10)/cdr2_mw,
         cdr3_percentextcoef = (cdr3_molextcoef*10)/cdr3_mw,
         fr1_percentextcoef = (fr1_molextcoef*10)/fr1_mw,
         fr2_percentextcoef = (fr2_molextcoef*10)/fr2_mw,
         fr3_percentextcoef = (fr3_molextcoef*10)/fr3_mw,
         fr4_percentextcoef = (fr4_molextcoef*10)/fr4_mw,
         AbChain_cysbridges_percenextcoef = (AbChain_cysbridges_molextcoef*10)/AbChain_mw,
         cdr1_cysbridges_percenextcoef = (cdr1_cysbridges_molextcoef*10)/cdr1_mw,
         cdr2_cysbridges_percenextcoef = (cdr2_cysbridges_molextcoef*10)/cdr2_mw,
         cdr3_cysbridges_percenextcoef = (cdr3_cysbridges_molextcoef*10)/cdr3_mw,
         fr1_cysbridges_percenextcoef = (fr1_cysbridges_molextcoef*10)/fr1_mw,
         fr2_cysbridges_percenextcoef = (fr2_cysbridges_molextcoef*10)/fr2_mw,
         fr3_cysbridges_percenextcoef = (fr3_cysbridges_molextcoef*10)/fr3_mw,
         fr4_cysbridges_percenextcoef = (fr4_cysbridges_molextcoef*10)/fr4_mw)

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

saveRDS(calculated_parameters, glue("{output_file}.RDS"))
fwrite(calculated_parameters, glue("{output_file}.csv"))
