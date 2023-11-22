library(dlookr)      # for exploratory data analysis and imputation
library(visdat)      # for visualizing NAs
library(plotly)      # for interactive visualization
library(missRanger)  # to generate and impute NAs
library(tidyverse)   
library(reshape2)    
library(ggrepel)



#loading and oganising predicition data (MLR from Jah, PLM from Matteo)----

PLM_predictions = read_csv("~/developability_project/plotting_dataframes/from_matteo/predictability/embs_linear_regressors_clipped.csv") %>%
  clean_names() %>% 
  pivot_longer(-c(x1), names_to = "sample_size", values_to = "rsquare") %>%
  rename(DP = x1) %>%
  mutate(sample_size = str_remove_all(sample_size, "x")) %>%
  mutate(rsquare = str_remove(rsquare, "\\[")) %>%
  mutate(rsquare = str_remove(rsquare, "\\]")) %>%
  mutate(rsquare = str_replace_all(rsquare,"  ", " ")) %>%
  mutate(rsquare = str_replace_all(rsquare," ", "_")) %>%
  mutate(rsquare = str_replace_all(rsquare,"_______|______|_____|____|___|__", "_")) %>%
  #filter(DP == "AbChain_mw") %>%
  #filter(sample_size == "100") %>% pull(rsquare)
  #mutate(rsquare = str_replace_all(rsquare,quant("_{2,8}"), "_")) %>%
  #mutate(rsquare = str_replace_all(rsquare,"__", "_"))  %>%
  separate(rsquare, as.character(c(1:20)), sep = "_") %>% 
  pivot_longer(-c(DP,sample_size), names_to = "replica", values_to = "rsquare_single") %>%
  mutate(rsquare_single = as.numeric(rsquare_single)) %>%
  mutate(level = case_when(str_detect(DP, "AbChain") ~ "sequence",
                           str_detect(DP, "AbStruc") ~ "structure"))



PLM_predictions_means = PLM_predictions %>%
  group_by(sample_size,DP) %>%
  summarise(mean_rsquare = mean(rsquare_single), sd_rsquare = sd(rsquare_single)) %>%
  mutate(level = case_when(str_detect(DP, "AbChain") ~ "sequence",
                           str_detect(DP, "AbStruc") ~ "structure")) %>% 
  #mutate(sample_size = as.double(sample_size)) %>%
  arrange(sample_size) %>%
  mutate(level = str_to_title(level))


PLM_predictions_means$sample_size = factor(PLM_predictions_means$sample_size, levels = c( "50","100","500","1000","10000","20000","100000"))
PLM_predictions_means_text = PLM_predictions_means %>%
  group_by(level,sample_size) %>%
  summarise(MD = median(mean_rsquare))


MLR_predictions = readRDS("~/developability_project/plotting_dataframes/from_jahn/predictability_june_2023/prediction_accuracies.rds")
    
MLR_predictions = MLR_predictions %>%
  rename(sample_size = n,
         rsquare_single = rsq,
         DP = parameter)



MLR_predictions_means = MLR_predictions %>%
  filter(dp_set_name == "comb_parameters_mwds") %>%
  mutate(rsquare_single = if_else(rsquare_single < 0, 0, rsquare_single)) %>% #crop the values to replace negative R squares with zero
  group_by(sample_size,DP) %>%
  summarise(mean_rsquare = mean(rsquare_single), sd_rsquare = sd(rsquare_single)) %>%
  mutate(level = case_when(str_detect(DP, "AbChain") ~ "sequence",
                           str_detect(DP, "AbStruc") ~ "structure")) %>% 
  #mutate(sample_size = as.double(sample_size)) %>%
  arrange(sample_size) %>%
  mutate(level = str_to_title(level)) 



y = MLR_predictions_means %>%
  group_by(level, sample_size) %>%
  summarise(mean_mean_rsquare = mean(mean_rsquare), sd_mean_rsquare = sd(mean_rsquare)) %>%
  mutate(method = "MLR") %>%
  mutate(sample_size = as.character(sample_size))

x_filtered = PLM_predictions_means %>%
  filter(DP %in% MLR_predictions_means$DP) %>% #ensure including same DPs in the comparison
  group_by(level, sample_size) %>%
  summarise(mean_mean_rsquare = mean(mean_rsquare), sd_mean_rsquare = sd(mean_rsquare)) %>%
  mutate(method = "PLM") 


joint_MLR_PLM_sample_size_filtered = rbind(x_filtered,y) %>%
  filter(sample_size != "100000") 

joint_MLR_PLM_sample_size_filtered = joint_MLR_PLM_sample_size_filtered %>%
  rename(Embeddings = method) %>%
  mutate(Embeddings = str_replace_all( Embeddings, pattern = "MLR", replacement = "DPLs")) 

#Figure 6B-----
ggplot(joint_MLR_PLM_sample_size_filtered, aes(group = Embeddings, x = sample_size, y = mean_mean_rsquare, color = Embeddings)) +
  facet_grid(cols = vars(level)) +
  geom_line()+
  geom_point(position=position_dodge(0.05))+
  geom_errorbar(aes(ymin = mean_mean_rsquare - sd_mean_rsquare, 
                    ymax = mean_mean_rsquare + sd_mean_rsquare),
                width=.2,
                position=position_dodge(0.05)) +
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"), strip.text = element_text(size = 12),
        legend.position = "top", legend.justification = "left",
        axis.text = element_text(size = 10)) +
  ylab(expression(Mean~R^2)) + xlab("Sample size (antibodies)") +
  scale_colour_manual(values = c("#1b9e77", "#d95f02")) +
  #geom_text(aes(label = round(mean_mean_rsquare, digits = 2)), show.legend = F) +
  geom_text_repel(aes(label = round(mean_mean_rsquare, digits = 2)), size = 2.5)

ggsave("~/developability_project/plots/predictability/DPs_PLM_seq_struc_line_comparison_corrected.pdf", width = 10, height = 5)
ggsave("~/developability_project/plots/predictability/DPs_PLM_seq_struc_line_comparison_corrected.png", width = 10, height = 5)


#generating MICRF predictions 





#Import MWDS sets------
AbStruc_mwds = readRDS('~/developability_project/plotting_dataframes/from_matteo/predictability/AbStruc_mwds.rds')
AbChain_mwds = readRDS('~/developability_project/plotting_dataframes/from_matteo/predictability/AbChain_mwds.rds')
AbChStruc_mwds = readRDS('~/developability_project/plotting_dataframes/from_matteo/predictability/AbChStruc_mwds.rds')



#remove immunology DPs and two structre DPs from the MWDS sets as to harmonise wiht Jahn and Matteo. 

AbChain_mwds_harmonised = AbChain_mwds[!AbChain_mwds %in% c("AbChain_num_strong_binders", "AbChain_min_rank", "AbChain_immunopeptide_regions_span")]
AbStruc_mwds_harmonised = AbStruc_mwds[!AbStruc_mwds %in% c("AbStruc_weak_hbonds", "AbStruc_pi_helices")]


# Import DP data------
setwd("~/developability_project/plotting_dataframes/from_matteo/predictability/reproducibility/")
AbChStruc_developability_normalised = np$load("human-heavy_metrics_normed.npy") %>%
  as_tibble()
# Import column names----
params = c(
  'AbChain_mw',
  'AbChain_length',
  'AbChain_av_resdiue_weight',
  'AbChain_aromatic_content',
  'AbChain_1_charge',
  'AbChain_2_charge',
  'AbChain_3_charge',
  'AbChain_4_charge',
  'AbChain_5_charge',
  'AbChain_6_charge',
  'AbChain_7_charge',
  'AbChain_8_charge',
  'AbChain_9_charge',
  'AbChain_10_charge',
  'AbChain_11_charge',
  'AbChain_12_charge',
  'AbChain_13_charge',
  'AbChain_14_charge',
  'AbChain_pI',
  'AbChain_molextcoef',
  'AbChain_cysbridges_molextcoef',
  'AbChain_percentextcoef',
  'AbChain_cysbridges_percenextcoef',
  'AbChain_hydrophobicity',
  'AbChain_hmom',
  'AbChain_solubility',
  'AbChain_instaindex',
  'AbChain_aliphindex',
  'AbChain_tiny_content',
  'AbChain_small_content',
  'AbChain_aliphatic_content',
  'AbChain_nonpolar_content',
  'AbChain_polar_content',
  'AbChain_basic_content',
  'AbChain_acidic_content',
  'AbStruc_loops',
  'AbStruc_amide_amide_interactions',
  'AbStruc_amide_ring_interactions',
  'AbStruc_beta_bridges',
  'AbStruc_carbon_pi_interactions',
  'AbStruc_cation_pi_interactions',
  'AbStruc_donor_pi_interactions',
  'AbStruc_beta_strands',
  'AbStruc_threeten_helices',
  'AbStruc_alpha_helices',
  'AbStruc_pi_helices', #should be removed
  'AbStruc_total_interactions',
  'AbStruc_metsulphur_pi_interactions',
  'AbStruc_beta_bends',
  'AbStruc_beta_turns',
  'AbStruc_aromatic_interactions',
  'AbStruc_atom_atom_sum',
  'AbStruc_atom_plane_sum',
  'AbStruc_carbonyl_interactions',
  'AbStruc_steric_clashes',
  'AbStruc_covbonds',
  'AbStruc_mean_interaction_distance',
  'AbStruc_hbonds',
  'AbStruc_hydrophobic_interactions',
  'AbStruc_ibonds',
  'AbStruc_phi_angle',
  'AbStruc_plane_plane_interactions',
  'AbStruc_pbonds',
  'AbStruc_proximal_interactions',
  'AbStruc_psi_angle',
  'AbStruc_vdw_interactions',
  'AbStruc_vdw_clashes',
  'AbStruc_weak_hbonds', #should be removed
  'AbStruc_weak_pbonds',
  'AbStruc_folding_energy',
  'AbStruc_unfolded_charge',
  'AbStruc_folded_charge',
  'AbStruc_sasa',
  'AbStruc_schains_sasa',
  'AbStruc_lauer_SAP',
  'AbStruc_lauer_DI',
  'AbStruc_free_cys',
  'AbStruc_cys_bridges',
  'AbStruc_pcharge_hetrgen',
  'AbStruc_ncharge_hetrgen',
  'AbStruc_plane_group_interactions',
  'AbStruc_folded_pI',
  'AbStruc_unfolded_pI'
)


# Add column names and index column-----
colnames(AbChStruc_developability_normalised) = params
AbChStruc_developability_normalised = AbChStruc_developability_normalised %>%
  mutate(ml_rowid = 1:nrow(.)) %>%
  select(ml_rowid, everything())
save(AbChStruc_developability_normalised, file = '~/developability_project/environments_DataObjects/AbChStruc_developability_normalised.RData')
load("~/developability_project/environments_DataObjects/AbChStruc_developability_normalised.RData")

AbChStruc_developability_normalised = AbChStruc_developability_normalised %>%
  select(-c(AbStruc_weak_hbonds, AbStruc_pi_helices))

sample_sizes = c('50', '100', '500', '1000', '10000', '20000')
batch_number = as.character(c(1:20))

#generate a rowid vectors that van be used whenb we are sampling and batching the scaled data
for (i in 1:length(sample_sizes)){
  for (j in 1:length(batch_number)){
    batch_size_rowids = AbChStruc_developability_normalised %>%
      sample_n(as.numeric(sample_sizes[i])) %>%
      pull(ml_rowid)
    print(paste0(batch_number[j], "_", sample_sizes[i]))
    save(batch_size_rowids, file = paste0("~/developability_project/environments_DataObjects/predictability_sets/", batch_number[j], "_", sample_sizes[i],".RData"))
  }
}

AbChain_mwds = AbChain_mwds

rsq <- function (observed, predicted) {
  1 - (sum((observed - predicted)^2) / sum((observed - mean(observed))^2))
}

save(rsq, file = "~/developability_project/environments_DataObjects/rsq.RDaata")


#run for sequence MWDS paramers predictions:
mice_predictions = tibble(sample_size = character(), batch = character(), level = character(), dp_set = character(), rsq = double(), data_loss = character())

for (i in 1:length(sample_sizes)){
  for (j in 1:length(batch_number)){
    load(paste0("~/developability_project/environments_DataObjects/predictability_sets/", batch_number[j], "_", sample_sizes[i],".RData"))
    #batch_size_rowids
    subset = AbChStruc_developability_normalised %>%
      filter(ml_rowid %in% batch_size_rowids) %>%
      select(-ml_rowid)
    #determine the number of values to be removed from the dataset (2%, then run for 4%)
    nmiss <- round(0.02*(ncol(subset)*nrow(subset))) 
    #nmiss <- round(0.04*(ncol(subset)*nrow(subset)))
    #original seq matrix
    subset_matrix_original_seq = subset %>%
      select(all_of(c(AbChain_mwds_harmonised))) %>% 
      as.matrix()
    #original struc matrix
    subset_matrix_original_struc = subset %>%
      select(all_of(c(AbStruc_mwds_harmonised))) %>% 
      as.matrix()
    subset_matrix_original =  cbind(subset_matrix_original_seq, subset_matrix_original_struc)
    #decide which parameters to remove (seq or struc)
    #remove seq DPs
    subset_matrix_missing = subset_matrix_original_seq
    positions_to_remove <- sample(nrow(subset_matrix_missing)*ncol(subset_matrix_missing), nmiss)
    subset_matrix_missing[positions_to_remove] <- NA
    
    observed_values <- subset_matrix_original_seq[positions_to_remove]
    
    #missing = sequence missing + original structure 
    subset_matrix_missing = cbind(subset_matrix_missing,subset_matrix_original_struc)
    dim(subset_matrix_missing)
    #it is working --> good!

    subset_matrix_missing_imputet = missRanger(
      as_tibble(subset_matrix_missing), 
      formula = . ~ . ,
      num.trees = 100, # The official page of the algorithm states that random forest does not overfit, and you can use as much trees as you want
      verbose = 1, 
      seed = 111,  
      returnOOB = F,
      num.threads = 10)
    
    subset_matrix_missing_imputet = subset_matrix_missing_imputet %>%
      as_tibble() %>%
      select(AbChain_mwds_harmonised) %>% #because we want the sequence parameters 
      as.matrix()
    
    predicted_values <- subset_matrix_missing_imputet[positions_to_remove]
    rsq_value = rsq(observed_values, predicted_values)
    
    mice_predictions = add_row(mice_predictions, sample_size = sample_sizes[i], batch = batch_number[j], 
                               level = "Sequence", dp_set = "comb_parameters_mwds", 
                               rsq = rsq_value, data_loss = "2%")
                               #rsq = rsq_value, data_loss = "4%")
    rm(subset_matrix_missing_imputet)
    rm(subset_matrix_missing)
    print(sample_sizes[i])
    print(batch_number[j])
  }
}

save(mice_predictions, file = "~/developability_project/environments_DataObjects/mice_predictions.RData")


#repeat the same for loop to do the structure DPs. 

for (i in 1:length(sample_sizes)){
  for (j in 1:length(batch_number)){
    load(paste0("~/developability_project/environments_DataObjects/predictability_sets/", batch_number[j], "_", sample_sizes[i],".RData"))
    #batch_size_rowids
    subset = AbChStruc_developability_normalised %>%
      filter(ml_rowid %in% batch_size_rowids) %>%
      select(-ml_rowid)
    #determine the number of values to be removed from the dataset (2%, then run for 4%)
    #nmiss <- round(0.02*(ncol(subset)*nrow(subset))) 
    nmiss <- round(0.04*(ncol(subset)*nrow(subset))) 
    #original seq matrix
    subset_matrix_original_seq = subset %>%
      select(all_of(c(AbChain_mwds_harmonised))) %>% 
      as.matrix()
    dim(subset_matrix_original_seq)
    #original struc matrix
    subset_matrix_original_struc = subset %>%
      select(all_of(c(AbStruc_mwds_harmonised))) %>% 
      as.matrix()
    dim(subset_matrix_original_struc)
    subset_matrix_original =  cbind(subset_matrix_original_seq, subset_matrix_original_struc)
    #dim(subset_matrix_original)
    #decide which parameters to remove (seq or struc)
    #remove seq DPs
    subset_matrix_missing = subset_matrix_original_struc
    positions_to_remove <- sample(nrow(subset_matrix_missing)*ncol(subset_matrix_missing), nmiss)
    subset_matrix_missing[positions_to_remove] <- NA
    
    observed_values <- subset_matrix_original_struc[positions_to_remove]
    
    #missing = sequence missing + original structure 
    subset_matrix_missing = cbind(subset_matrix_original_seq,subset_matrix_missing)
    #dim(subset_matrix_missing)
    #it is working --> good!
    #sum(is.na(subset_matrix_original))
    #sum(is.na(subset_matrix_missing))
    
    subset_matrix_missing_imputet = missRanger(
      as_tibble(subset_matrix_missing), 
      formula = . ~ . ,
      num.trees = 100, # The official page of the algorithm states that random forest does not overfit, and you can use as much trees as you want
      verbose = 1, 
      seed = 111,  
      returnOOB = F,
      num.threads = 10)
    
    subset_matrix_missing_imputet = subset_matrix_missing_imputet %>%
      as_tibble() %>%
      select(all_of(AbStruc_mwds_harmonised)) %>% #because we want the struc parameters 
      as.matrix()
    
    predicted_values <- subset_matrix_missing_imputet[positions_to_remove]
    rsq_value = rsq(observed_values, predicted_values)
    
    mice_predictions = add_row(mice_predictions, sample_size = sample_sizes[i], batch = batch_number[j], 
                               level = "Structure", dp_set = "comb_parameters_mwds", 
                               #rsq = rsq_value, data_loss = "2%")
                               rsq = rsq_value, data_loss = "4%")
    rm(subset_matrix_missing_imputet)
    rm(subset_matrix_missing)
    
    print(sample_sizes[i])
    print(batch_number[j])
  }
}

#run for 2%, repeat for 4%
save(mice_predictions, file = "~/developability_project/environments_DataObjects/mice_predictions.RData")


#then repeat for NON-MWDS DPs

#make any prediction below 0 equal to zero 

mice_predictions = mice_predictions %>%
  mutate(rsq = case_when(rsq < 0 ~ 0,TRUE ~ rsq)) #%>%
#View()

mice_predictions_mean = mice_predictions %>%
  group_by(sample_size, level, data_loss, dp_set) %>%
  summarise(mean_rsq = mean(rsq), sd_rsq = sd(rsq))

mice_predictions_mean$sample_size = factor(mice_predictions_mean$sample_size, levels = c( "50","100","500","1000","10000","20000"))

#Figure 6C-----

ggplot(mice_predictions_mean %>% filter(dp_set == "comb_parameters_mwds"), 
       aes(group = data_loss, x = sample_size, y = mean_rsq, color = data_loss)) +
  facet_grid(cols = vars(level)) +
  geom_line()+
  geom_point(position=position_dodge(0.05))+
  geom_errorbar(aes(ymin = mean_rsq - sd_rsq, 
                    ymax = mean_rsq + sd_rsq),
                width=.2,
                position=position_dodge(0.05)) +
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"), strip.text = element_text(size = 12),
        legend.position = "top", legend.justification = "left",
        axis.text = element_text(size = 10)) +
  ylab(expression(Mean~R^2)) + xlab("Training set size (antibodies)") +
  scale_colour_manual(values = c("deepskyblue", "dodgerblue4"),
                      guide_legend(title="Data loss")) +
  geom_text_repel(aes(label = round(mean_rsq, digits = 2)), size = 2.5, show.legend = F) +
  ylim(c(NA,1))


ggsave("~/Documents/Work/Oslo/Antibody_developability/plots/predictability/twopercent_fourpercent_line_comparison_MICE.pdf", width = 10, height = 5)
ggsave("~/Documents/Work/Oslo/Antibody_developability/plots/predictability/twopercent_fourpercent_line_comparison_MICE.png", width = 10, height = 5)







#Supp Figure 14A top panel-----

MLR_predictions_means_set_comparisons = MLR_predictions %>%
  #filter(dp_set_name == "comb_parameters_mwds") %>%
  mutate(rsquare_single = if_else(rsquare_single < 0, 0, rsquare_single)) %>% #crop the values to replace negative R squares with zero
  group_by(sample_size,dp_set_name,DP) %>%
  summarise(mean_rsquare = mean(rsquare_single), sd_rsquare = sd(rsquare_single)) %>%
  mutate(level = case_when(str_detect(DP, "AbChain") ~ "sequence",
                           str_detect(DP, "AbStruc") ~ "structure")) %>% 
  #mutate(sample_size = as.double(sample_size)) %>%
  #arrange(sample_size) %>%
  mutate(level = str_to_title(level)) %>%
  filter(dp_set_name %in% c("comb_parameters_mwds","comb_parameters_subsample")) %>%
  group_by(level, dp_set_name, sample_size) %>%
  summarise(mean_mean_rsquare = mean(mean_rsquare), sd_mean_rsquare = sd(mean_rsquare))


MLR_predictions_means_set_comparisons = MLR_predictions_means_set_comparisons %>%
  mutate(sample_size = as.character(sample_size))

MLR_predictions_means_set_comparisons$sample_size = factor(MLR_predictions_means_set_comparisons$sample_size, levels = c( "50","100","500","1000","10000","20000"))

MLR_predictions_means_set_comparisons = MLR_predictions_means_set_comparisons %>%
  mutate(dp_set_name = str_replace_all(dp_set_name, "comb_parameters_mwds", "MWDS")) %>%
  mutate(dp_set_name = str_replace_all(dp_set_name, "comb_parameters_subsample", "non-MWDS")) 


ggplot(MLR_predictions_means_set_comparisons, aes(group = dp_set_name, x = sample_size, y = mean_mean_rsquare, color = dp_set_name)) +
  facet_grid(cols = vars(level)) +
  geom_line()+
  geom_point(position=position_dodge(0.05))+
  geom_errorbar(aes(ymin = mean_mean_rsquare - sd_mean_rsquare, 
                    ymax = mean_mean_rsquare + sd_mean_rsquare),
                width=.2,
                position=position_dodge(0.05)) +
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"), strip.text = element_text(size = 12),
        legend.position = "top", legend.justification = "left",
        axis.text = element_text(size = 10)) +
  ylab(expression(Mean~R^2)) + xlab("Sample size (antibodies)") +
  scale_colour_manual(values = c("#d8b365", "#5ab4ac"),
                      guide_legend(title="Dominance state")) +
  geom_text_repel(aes(label = round(mean_mean_rsquare, digits = 2)), size = 2.5)


ggsave("~/developability_project/plots/predictability/DPs_MWDS_nonMWDS.pdf", width = 10, height = 4)
ggsave("~/developability_project/plots/predictability/DPs_MWDS_nonMWDS.png", width = 10, height = 4)


#Supp Figure 14A bottom panel-----

xx_filtered = PLM_predictions_means %>%
  filter(DP %in% non_MWDS_DPs) %>% #ensure including same DPs in the comparison
  group_by(level, sample_size) %>%
  summarise(mean_mean_rsquare = mean(mean_rsquare), sd_mean_rsquare = sd(mean_rsquare)) %>%
  mutate(method = "non-MWDS") 

rbind(x_filtered, xx_filtered) %>%
  filter(sample_size != "100000") %>%
  ggplot(aes(group = method, x = sample_size, y = mean_mean_rsquare, color = method)) +
  facet_grid(cols = vars(level)) +
  geom_line()+
  geom_point(position=position_dodge(0.05))+
  geom_errorbar(aes(ymin = mean_mean_rsquare - sd_mean_rsquare, 
                    ymax = mean_mean_rsquare + sd_mean_rsquare),
                width=.2,
                position=position_dodge(0.05)) +
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"), strip.text = element_text(size = 12),
        legend.position = "top", legend.justification = "left",
        axis.text = element_text(size = 10)) +
  ylab(expression(Mean~R^2)) + xlab("Sample size (antibodies)") +
  scale_colour_manual(values = c("#d8b365", "#5ab4ac"),
                      guide_legend(title="Dominance state")) +
  geom_text_repel(aes(label = round(mean_mean_rsquare, digits = 2)), size = 2.5)

ggsave("~/developability_project/plots/predictability/PLM_MWDS_nonMWDS.pdf", width = 10, height = 4)
ggsave("~/developability_project/plots/predictability/PLM_MWDS_nonMWDS.png", width = 10, height = 4)

#Supp Figure 14B------
ggplot(mice_predictions_mean %>% filter(dp_set == "comb_parameters_subsample"), 
       aes(group = data_loss, x = sample_size, y = mean_rsq, color = data_loss)) +
  facet_grid(cols = vars(level)) +
  geom_line()+
  geom_point(position=position_dodge(0.05))+
  geom_errorbar(aes(ymin = mean_rsq - sd_rsq, 
                    ymax = mean_rsq + sd_rsq),
                width=.2,
                position=position_dodge(0.05)) +
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"), strip.text = element_text(size = 12),
        legend.position = "top", legend.justification = "left",
        axis.text = element_text(size = 10)) +
  ylab(expression(Mean~R^2)) + xlab("Training set size (antibodies)") +
  scale_colour_manual(values = c("deepskyblue", "dodgerblue4"),
                      guide_legend(title="Data loss")) +
  geom_text_repel(aes(label = round(mean_rsq, digits = 2)), size = 2.5, show.legend = F) +
  ylim(c(NA,1))



ggsave("~/Documents/Work/Oslo/Antibody_developability/plots/predictability/twopercent_fourpercent_line_comparison_MICE_non_MWDS.pdf", width = 10, height = 5)
ggsave("~/Documents/Work/Oslo/Antibody_developability/plots/predictability/twopercent_fourpercent_line_comparison_MICE_non_MWDS.png", width = 10, height = 5)

