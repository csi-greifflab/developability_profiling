#Author: Habib Bashour

library(reshape2)
library(tidyverse)
library(tidyverse)
library(reticulate)
#Figure 7A----- (correlation distance dendrograms)-----
#load correlation matrices generated as previously explained in the redundancy section for the native and human-engineered antibodies

loadRData = function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

setwd("~/developability_project/environments_DataObjects/parameter_correlation_matrices_pearson/fv/")
#run alternatively (once for sequence DPs and once for structure DPs)

list_RData = dir(pattern = "AbChain")
#list_RData = dir(pattern = "AbStruc")


#get one dataframe for native - sequence only

#run the loop for i = 1 until drop_na() then run the following
ordered_cor_df = sample_pairwise_data %>% select(-value)

i = 1
for(i in 1:length(list_RData)){
  matrix = loadRData(paste0("~/developability_project/environments_DataObjects/parameter_correlation_matrices_pearson/fv/",list_RData[i]))
  file_name = str_remove_all(list_RData[i], "AbChain_parameter_cormatrix_fv_|.RData")
  #file_name = str_remove_all(list_RData[i], "AbStruc_parameter_cormatrix_fv_|.RData")
  diag(matrix) = NA
  matrix[lower.tri(matrix)] = NA
  sample_pairwise_data = matrix %>% 
    melt() %>% 
    as_tibble() %>% 
    drop_na() 
  
  #rename() %>%
  #unite(pair,c("Var1", "Var2"), sep = "_")
  sample_pairwise_data = setNames(sample_pairwise_data, c("Var1", "Var2",file_name))
  ordered_cor_df = full_join(ordered_cor_df,sample_pairwise_data) #full_join so we do not lose columns when a matrix is missing the pair
}
AbChain_ordered_cor_df = ordered_cor_df
#AbStruc_ordered_cor_df = ordered_cor_df


#use the same running instructions as above for mAbs
setwd("~/developability_project/environments_DataObjects/parameter_correlation_matrices_pearson/mAbs/")
list_RData = dir(pattern = "AbChain")
#list_RData = dir(pattern = "AbStruc")
ordered_cor_df_mAbs = sample_pairwise_data %>% select(-value)

i = 1
for(i in 1:length(list_RData)){
  matrix = loadRData(paste0("~/developability_project/environments_DataObjects/parameter_correlation_matrices_pearson/mAbs/",list_RData[i]))
  #file_name = str_remove_all(list_RData[i], "AbChain_parameter_cormatrix_fv_|.RData")
  file_name = str_remove_all(list_RData[i], "AbStruc_parameter_cormatrix_fv_|.RData")
  diag(matrix) = NA
  #matrix[lower.tri(matrix)] = NA
  sample_pairwise_data = matrix %>% 
    melt() %>% 
    as_tibble() %>% 
    mutate(dp_paires = paste0(Var1, "_", Var2)) %>%
    filter(dp_paires %in% dp_pairs_order_AbChain) %>%
    #filter(dp_paires %in% dp_pairs_order_AbStruc) %>%
    select(-dp_paires) %>%
    drop_na()
  
  sample_pairwise_data = setNames(sample_pairwise_data, c("Var1", "Var2",file_name))
  ordered_cor_df_mAbs = left_join(ordered_cor_df_mAbs,sample_pairwise_data) #full_join so we do not lose columns when a matrix is missing the pair
}

AbChain_ordered_cor_df_mAbs = ordered_cor_df_mAbs
#AbStruc_ordered_cor_df_mAbs = ordered_cor_df_mAbs



#now get for PADs

setwd("~/developability_project/environments_DataObjects/parameter_correlation_matrices_pearson/PAD/")
list_RData = dir(pattern = "AbChain")
#list_RData = dir(pattern = "AbStruc")


#run the loop for i = 1 until drop_na() then run the following
ordered_cor_df_PAD = sample_pairwise_data %>% select(-value)

i = 1
for(i in 1:length(list_RData)){
  matrix = loadRData(paste0("~/developability_project/environments_DataObjects/parameter_correlation_matrices_pearson/PAD/",list_RData[i]))
  file_name = str_remove_all(list_RData[i], "AbChain_parameter_cormatrix_fv_|.RData")
  #file_name = str_remove_all(list_RData[i], "AbStruc_parameter_cormatrix_fv_|.RData")
  diag(matrix) = NA
  sample_pairwise_data = matrix %>% 
    melt() %>% 
    as_tibble() %>%
    mutate(dp_paires = paste0(Var1, "_", Var2)) %>%
    filter(dp_paires %in% dp_pairs_order_AbChain) %>%
    #filter(dp_paires %in% dp_pairs_order_AbStruc) %>%
    select(-dp_paires) %>%
    drop_na()
  sample_pairwise_data = setNames(sample_pairwise_data, c("Var1", "Var2",file_name))
  ordered_cor_df_PAD = full_join(ordered_cor_df_PAD,sample_pairwise_data) #full_join so we do not lose columns when a matrix is missing the pair
}

AbChain_ordered_cor_df_PAD = ordered_cor_df_PAD
#AbStruc_ordered_cor_df_PAD  = ordered_cor_df_PAD


#now get for hum_mice

setwd("~/developability_project/environments_DataObjects/parameter_correlation_matrices_pearson/humanisd_mice/")
#run the loop for i = 1 until drop_na() then run the following
matrix = loadRData("~/developability_project/environments_DataObjects/parameter_correlation_matrices_pearson/humanisd_mice/AbChain_parameter_cormatrix_fv_humice.RData")
#matrix = loadRData("~/developability_project/environments_DataObjects/parameter_correlation_matrices_pearson/humanisd_mice/AbStruc_parameter_cormatrix_fv_humice.RData")
diag(matrix) = NA
sample_pairwise_data = matrix %>% 
  melt() %>% 
  as_tibble() %>% 
  mutate(dp_paires = paste0(Var1, "_", Var2)) %>%
  filter(dp_paires %in% dp_pairs_order_AbChain) %>%
  #filter(dp_paires %in% dp_pairs_order_AbStruc) %>%
  select(-dp_paires) 
drop_na()

AbChain_ordered_cor_df_humice = setNames(sample_pairwise_data, c("Var1", "Var2","kymouse"))
AbStruc_ordered_cor_df_humice = setNames(sample_pairwise_data, c("Var1", "Var2","kymouse"))



AbChain_all_ordered_cor_df = left_join(AbChain_ordered_cor_df, AbChain_ordered_cor_df_PAD) %>%
  full_join(AbChain_ordered_cor_df_mAbs) %>% 
  full_join(AbChain_ordered_cor_df_humice) %>% 
  drop_na() %>%
  select(-pair)



#calculate correlation distance 

scale_manual_fun_dist_optimised = function(x){
  (x - mean(x)) / sd(x)
}


my_get_dist = function(x, method = "euclidean", stand = FALSE){
  if (stand) 
    x <- x %>% mutate_at(scale_manual_fun_dist_optimised, .vars = vars(-c(Var1,Var2))) #this is required 
  x <- x %>% select(-c(Var1,Var2))
  if (method %in% c("pearson", "spearman", "kendall")) {
    res.cor <- stats::cor(x, method = method, use = "pairwise.complete.obs") #I don't want to flip
    res.dist <- stats::as.dist(1 - res.cor)
  }
  else res.dist <- stats::dist(x, method = method, ...)
  res.dist
  
}


AbChain_all_dist = my_get_dist(AbChain_all_ordered_cor_df, method = "pearson", stand = TRUE) 
AbChain_all_dend <- AbChain_all_dist %>% hclust %>% as.dendrogram

#Figure 7A top panel-
pdf("~/developability_project/plots/cor_correlation/dendrograms/pearson_distance_pairwise_all_AbChain_uncoloured.pdf", width = 9, height = 8)
par(cex=1)
AbChain_all_dend %>%  
  set("branches_lwd", 2) %>%
  plot(lab = "correlation distance") 
dev.off()




AbStruc_all_ordered_cor_df = left_join(AbStruc_ordered_cor_df, AbStruc_ordered_cor_df_PAD) %>%
  full_join(AbStruc_ordered_cor_df_mAbs) %>% 
  full_join(AbStruc_ordered_cor_df_humice) %>% 
  drop_na()


AbStruc_all_dist = my_get_dist(AbStruc_all_ordered_cor_df, method = "pearson", stand = TRUE) 
AbStruc_all_dend <- AbStruc_all_dist %>% hclust %>% as.dendrogram

#Figure 7A bottom panel-
pdf("~/developability_project/plots/cor_correlation/dendrograms/pearson_distance_pairwise_all_AbStru_uncoloured.pdf", width = 9, height = 8)
par(cex=1)
AbStruc_all_dend %>%  
  set("branches_lwd", 2) %>%
  plot(ylab = "correlation distance") 
dev.off()



#Figure 7B----(overlaying datasets in the developability space)-----

#the heavy chains

np <- import("numpy")

ed_ditance_matrix_IgM <- np$load("sample_embs_dists.npy")
#directory of heavy human
setwd("~/developability_project/plotting_dataframes/from_matteo/Native_non-Native_Comparison_complete/")
list_PLM_proj = dir(pattern = "*esm-1v-1_proj.npy")
list_DPL_proj = dir(pattern = "*mwds_proj.npy")
list_metadata = dir(pattern = "*_metadata.csv")



#run for i = 1 to create the empty tibble

DPL_PLM_PCAs_full = cbind(metadata,PLM, DPL) %>%
  as_tibble() %>%
  filter(dataset == "habib")


for(i in 1:length(list_metadata)){
  setwd("~/developability_project/plotting_dataframes/from_matteo/Native_non-Native_Comparison_complete/")
  dataset = list_metadata[i] %>% str_remove("v3_heavy-human-developability-")  %>% str_remove("_metadata.csv")
  metadata = read_csv(paste0("~/developability_project/plotting_dataframes/from_matteo/Native_non-Native_Comparison_complete/", list_metadata[i])) %>% 
    select(rowid, starts_with("AbOrig"))  %>%
    mutate(dataset = dataset)
  PLM =  np$load(list_PLM_proj[i]) %>%
    as_tibble() %>%
    select(V1,V2) %>%
    rename(PC1_PLM = V1, PC2_PLM = V2) 
  
  DPL =  np$load(list_DPL_proj[i]) %>%
    as_tibble() %>%
    rename(PC1_DPL = V1, PC2_DPL = V2)
  
  single = cbind(metadata,PLM, DPL) %>%
    as_tibble()
  DPL_PLM_PCAs_full = rbind(single, DPL_PLM_PCAs_full)
  rm(DPL,PLM, metadata, single)
}

save(DPL_PLM_PCAs_full, file = "~/developability_project/environments_DataObjects/DPL_PLM_PCAs_full.RData")
load("~/developability_project/environments_DataObjects/DPL_PLM_PCAs_full.RData")



####the following three plots form the top three panels of Figure 7B

ggplot(DPL_PLM_PCAs_full %>% filter(dataset == "native"), aes(x = PC1_DPL, y = PC2_DPL))+
  stat_binhex(aes(x = PC1_DPL, y = PC2_DPL), bins = 100)+ 
  geom_point(data = DPL_PLM_PCAs_full %>% filter(dataset == "thera"),
             aes(x = PC1_DPL, y = PC2_DPL), 
             shape = 21,
             fill = "#5e3c99",
             colour = "white",
             size = 2, alpha = 0.7) +
  theme_bw()+
  theme(#legend.position = "top", legend.justification = "left",legend.key.width= unit(1.5, 'cm'), [plot once with this line to extract the legend the legend]
    legend.position = "none",
    text = element_text(family="Helvetica"),
    axis.title = element_text(size = 12)) +
  labs(fill = "Native antibody count") +
  xlab("PC1 (9%)") + ylab("PC2 (8%)") #Source: @Matteo Pariset

ggsave("~/developability_project/plots/PCAs/integration_PCA_DPL_embeddings_native_and_mabs_hex.png", height = 6, width = 6)
ggsave("~/developability_project/plots/PCAs/integration_PCA_DPL_embeddings_native_and_mabs_hex.pdf", height = 6, width = 6)

p = ggplot(DPL_PLM_PCAs_full %>% filter(dataset == "native"), aes(x = PC1_DPL, y = PC2_DPL))+
  stat_binhex(aes(x = PC1_DPL, y = PC2_DPL), bins = 100)+ 
  geom_point(data = DPL_PLM_PCAs_full %>% filter(dataset == "humanised"),
             aes(x = PC1_DPL, y = PC2_DPL), colour = "#f5c37c", 
             size = 0.5, alpha = 0.5) +
  theme_bw()+
  theme(#legend.position = "top", legend.justification = "left",legend.key.width= unit(1.5, 'cm'),
    legend.position = "none",
    text = element_text(family="Helvetica"),
    axis.title = element_text(size = 12)) +
  labs(fill = "Native antibody count") +
  xlab("PC1 (9%)") + ylab("PC2 (8%)")

ggsave(p, file = "~/developability_project/plots/PCAs/integration_PCA_DPL_embeddings_native_and_kymouse_hex.png", height = 6, width = 6)
ggsave(p, file = "~/developability_project/plots/PCAs/integration_PCA_DPL_embeddings_native_and_kymouse_hex.pdf", height = 6, width = 6)


p = ggplot(DPL_PLM_PCAs_full %>% filter(dataset == "native"), aes(x = PC1_DPL, y = PC2_DPL))+
  stat_binhex(aes(x = PC1_DPL, y = PC2_DPL), bins = 100)+ 
  geom_point(data = DPL_PLM_PCAs_full %>% filter(dataset == "patented"),
             aes(x = PC1_DPL, y = PC2_DPL), colour = "#b2abd2", 
             size = 0.5, alpha = 0.5) +
  theme_bw()+
  theme(#legend.position = "top", legend.justification = "left",legend.key.width= unit(1.5, 'cm'),
    legend.position = "none",
    text = element_text(family="Helvetica"),
    axis.title = element_text(size = 12)) +
  labs(fill = "Native antibody count") +
  xlab("PC1 (9%)") + ylab("PC2 (8%)")


ggsave(p, file = "~/developability_project/plots/PCAs/integration_PCA_DPL_embeddings_native_and_pads_hex.png", height = 6, width = 6)
ggsave(p, file = "~/developability_project/plots/PCAs/integration_PCA_DPL_embeddings_native_and_pads_hex.pdf", height = 6, width = 6)


#the light chains 

setwd("~/developability_project/plotting_dataframes/from_matteo/native_non-native_comparison_light/")
list_PLM_proj = dir(pattern = "*esm-1v-1_proj.npy") #no himanised mouse --> three files
list_DPL_proj = dir(pattern = "*mwds_proj.npy") 
list_metadata = dir(pattern = "*_metadata.csv")


for(i in 1:length(list_metadata)){
  setwd("~/developability_project/plotting_dataframes/from_matteo/native_non-native_comparison_light/")
  dataset = list_metadata[i] %>% str_remove("v3_light-human-developability-")  %>% str_remove("_metadata.csv")
  metadata = read_csv(paste0("~/developability_project/plotting_dataframes/from_matteo/native_non-native_comparison_light/", list_metadata[i])) %>% 
    select(rowid, starts_with("AbOrig"))  %>%
    mutate(dataset = dataset)
  PLM =  np$load(list_PLM_proj[i]) %>%
    as_tibble() %>%
    select(V1,V2) %>%
    rename(PC1_PLM = V1, PC2_PLM = V2) 
  
  DPL =  np$load(list_DPL_proj[i]) %>%
    as_tibble() %>%
    rename(PC1_DPL = V1, PC2_DPL = V2)
  
  #Matteo said that the columns are all in the same order:
  single = cbind(metadata,PLM, DPL) %>%
    as_tibble()
  DPL_PLM_PCAs_full_light = rbind(single, DPL_PLM_PCAs_full_light)
  rm(DPL,PLM, metadata, single)
}

#plotting 


save(DPL_PLM_PCAs_full_light, file = "~/developability_project/environments_DataObjects/DPL_PLM_PCAs_full_light.RData")


#the following two plots form the bottom panels of Figure 7B

ggplot(DPL_PLM_PCAs_full_light %>% filter(dataset == "native"), aes(x = PC1_DPL, y = PC2_DPL))+
  stat_binhex(aes(x = PC1_DPL, y = PC2_DPL), bins = 100)+ 
  geom_point(data = DPL_PLM_PCAs_full_light %>% filter(dataset == "thera"),
             aes(x = PC1_DPL, y = PC2_DPL), 
             shape = 21,
             fill = "#5e3c99",
             colour = "white",
             size = 2, alpha = 0.7) +
  theme_bw()+
  theme(#legend.position = "top", legend.justification = "left",legend.key.width= unit(1.5, 'cm'),
    legend.position = "none",
    text = element_text(family="Helvetica"),
    axis.title = element_text(size = 12)) +
  labs(fill = "Native antibody count") +
  xlab("PC1 (20%)") + ylab("PC2 (13%)")
?element_text

ggsave("~/developability_project/plots/PCAs/integration_PCA_DPL_embeddings_native_and_mabs_hex_legend.pdf", height = 6, width = 6)
ggsave("~/developability_project/plots/PCAs/integration_PCA_DPL_embeddings_native_and_mabs_hex_light.png", height = 6, width = 9)
ggsave("~/developability_project/plots/PCAs/integration_PCA_DPL_embeddings_native_and_mabs_hex_light.pdf", height = 6, width = 9)



p = ggplot(DPL_PLM_PCAs_full_light %>% filter(dataset == "native"), aes(x = PC1_DPL, y = PC2_DPL))+
  stat_binhex(aes(x = PC1_DPL, y = PC2_DPL), bins = 100)+ 
  geom_point(data = DPL_PLM_PCAs_full_light %>% filter(dataset == "patented"),
             aes(x = PC1_DPL, y = PC2_DPL), colour = "#b2abd2", 
             size = 0.5, alpha = 0.5) +
  theme_bw()+
  theme(#legend.position = "top", legend.justification = "left",legend.key.width= unit(1.5, 'cm'),
    legend.position = "none",
    text = element_text(family="Helvetica"),
    axis.title = element_text(size = 12)) +
  labs(fill = "Native antibody count") +
  xlab("PC1 (20%)") + ylab("PC2 (13%)")


#ggsave("~/developability_project/plots/PCAs/integration_PCA_DPL_embeddings_native_and_pads.png", height = 5, width = 5)
ggsave(p, file = "~/developability_project/plots/PCAs/integration_PCA_DPL_embeddings_native_and_pads_hex_light.png", height = 6, width = 9)
ggsave(p, file = "~/developability_project/plots/PCAs/integration_PCA_DPL_embeddings_native_and_pads_hex_light.pdf", height = 6, width = 9)




#repeat the same plots but this time for PLM embeddings (for supp figures)----------



ggplot(DPL_PLM_PCAs_full %>% filter(dataset == "native"), aes(x = PC1_PLM, y = PC2_PLM))+
  stat_binhex(aes(x = PC1_PLM, y = PC2_PLM), bins = 100)+ 
  geom_point(data = DPL_PLM_PCAs_full %>% filter(dataset == "thera"),
             aes(x = PC1_PLM, y = PC2_PLM), 
             shape = 21,
             fill = "#5e3c99",
             colour = "white",
             size = 2, alpha = 0.7) +
  theme_bw()+
  theme(#legend.position = "top", legend.justification = "left",legend.key.width= unit(1.5, 'cm'), #run once for the legend
    legend.position = "none",
    text = element_text(family="Helvetica"),
    axis.title = element_text(size = 12)) +
  labs(fill = "Native antibody count") +
  xlab("PC1 (39%)") + ylab("PC2 (13%)")

ggsave("~/developability_project/plots/PCAs/integration_PCA_PLM_embeddings_native_and_mabs_hex.png", height = 6, width = 6)
ggsave("~/developability_project/plots/PCAs/integration_PCA_PLM_embeddings_native_and_mabs_hex.pdf", height = 6, width = 6)

p = ggplot(DPL_PLM_PCAs_full %>% filter(dataset == "native"), aes(x = PC1_PLM, y = PC2_PLM))+
  stat_binhex(aes(x = PC1_PLM, y = PC2_PLM), bins = 100)+ 
  geom_point(data = DPL_PLM_PCAs_full %>% filter(dataset == "humanised"),
             aes(x = PC1_PLM, y = PC2_PLM), colour = "#f5c37c", 
             size = 0.5, alpha = 0.5) +
  theme_bw()+
  theme(#legend.position = "top", legend.justification = "left",legend.key.width= unit(1.5, 'cm'),
    legend.position = "none",
    text = element_text(family="Helvetica"),
    axis.title = element_text(size = 12)) +
  labs(fill = "Native antibody count") +
  xlab("PC1 (39%)") + ylab("PC2 (13%)")

#ggsave(p, file = "~/developability_project/plots/PCAs/integration_PCA_DPL_embeddings_native_and_kymouse.png", height = 5, width = 5)
ggsave(p, file = "~/developability_project/plots/PCAs/integration_PCA_PLM_embeddings_native_and_kymouse_hex.png", height = 6, width = 6)
ggsave(p, file = "~/developability_project/plots/PCAs/integration_PCA_PLM_embeddings_native_and_kymouse_hex.pdf", height = 6, width = 6)


p = ggplot(DPL_PLM_PCAs_full %>% filter(dataset == "native"), aes(x = PC1_PLM, y = PC2_PLM))+
  #geom_hex(aes(x = PC1_DPL, y = PC2_DPL), bins = 30, colour = "white")+ 
  #stat_density_2d(aes(x = PC1_DPL, y = PC2_DPL, fill  = after_stat(level)),geom = "polygon", colour= NA) +
  stat_binhex(aes(x = PC1_PLM, y = PC2_PLM), bins = 100)+ 
  geom_point(data = DPL_PLM_PCAs_full %>% filter(dataset == "patented"),
             aes(x = PC1_PLM, y = PC2_PLM), colour = "#b2abd2", 
             size = 0.5, alpha = 0.5) +
  theme_bw()+
  theme(#legend.position = "top", legend.justification = "left",legend.key.width= unit(1.5, 'cm'),
    legend.position = "none",
    text = element_text(family="Helvetica"),
    axis.title = element_text(size = 12)) +
  labs(fill = "Native antibody count") +
  xlab("PC1 (39%)") + ylab("PC2 (13%)")


ggsave(p, file = "~/developability_project/plots/PCAs/integration_PCA_PLM_embeddings_native_and_pads_hex.png", height = 6, width = 6)
ggsave(p, file = "~/developability_project/plots/PCAs/integration_PCA_PLM_embeddings_native_and_pads_hex.pdf", height = 6, width = 6)


#check the role of isotype and IGHV in positioning in DPL space


DPL_PLM_PCAs_full$AbOrig_isotype = factor(DPL_PLM_PCAs_full$AbOrig_isotype, levels = c("IgD", "IgM", "IgG" ,"IgA" ,"IgE",
                                                                                       "IgH",
                                                                                       "IgG1","IgG2","IgG3","IgG4", "IgG4-G1"))  

nrow = DPL_PLM_PCAs_full %>% filter(dataset == "native") %>% nrow()
p = ggplot(DPL_PLM_PCAs_full %>% filter(dataset == "native") %>% sample_n(nrow), aes(x = PC1_DPL, y = PC2_DPL)) +
  geom_point(aes(x = PC1_DPL, y = PC2_DPL, colour = AbOrig_isotype), 
             size = 0.25, alpha = 1) +
  theme_bw()+
  xlab("PC1 (9%)") + ylab("PC2 (8%)") +
  scale_colour_manual(values = c("#7a0000", "#a80000", "#c10600", "#c45e05", "#c78219"),
                      name = "Isotype") +
  theme(axis.title = element_text(size = 12), legend.position = "top", legend.justification = "left",
        legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        text = element_text(family="Helvetica")) +
  guides(colour = guide_legend(override.aes = list(size=6)))

ggsave(p, file = "~/developability_project/plots/PCAs/integration_supporting_figure_by_isotype.png", height = 5.5, width = 5)
ggsave(p, file = "~/developability_project/plots/PCAs/integration_supporting_figure_by_isotype.pdf", height = 5.5, width = 5)

DPL_PLM_PCAs_full %>% filter(dataset == "native") %>% select(AbOrig_v_gene_prefix) %>% drop_na() #852,326/854,408 --> 2K has no v-gene annotations?

p = ggplot(DPL_PLM_PCAs_full %>% filter(dataset == "native") %>% sample_n(nrow), aes(x = PC1_DPL, y = PC2_DPL)) +
  geom_point(aes(x = PC1_DPL, y = PC2_DPL, colour = AbOrig_v_gene_prefix), 
             size = 0.35, alpha = 1) +
  theme_bw()+
  xlab("PC1 (9%)") + ylab("PC2 (8%)") +
  scale_colour_brewer(palette = "RdYlBu", name = "IGHV family") +
  theme(axis.title = element_text(size = 12), legend.position = "top", legend.justification = "left",
        legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        text = element_text(family="Helvetica")) +
  guides(#colour = "none",
    colour = guide_legend(override.aes = list(size=6)))

ggsave(p, file = "~/developability_project/plots/PCAs/integration_supporting_figure_by_vgene.png", height = 5.5, width = 5)

#check the role of isotype and IGHV in positioning in PLM space


p = ggplot(DPL_PLM_PCAs_full %>% filter(dataset == "native") %>% sample_n(nrow), aes(x = PC1_PLM, y = PC2_PLM)) +
  geom_point(aes(x = PC1_PLM, y = PC2_PLM, colour = AbOrig_v_gene_prefix), 
             size = 0.35, alpha = 1) +
  theme_bw()+
  xlab("PC1 (39%)") + ylab("PC2 (13%)") + 
  scale_colour_brewer(palette = "RdYlBu", name = "IGHV family") +
  theme(axis.title = element_text(size = 12), legend.position = "top", legend.justification = "left",
        legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        text = element_text(family="Helvetica")) +
  guides(#colour = "none",
    colour = guide_legend(override.aes = list(size=6)))

ggsave(p, file = "~/developability_project/plots/PCAs/integration_supporting_figure_PLM_by_vgene.png", height = 5.5, width = 5)

p = ggplot(DPL_PLM_PCAs_full %>% filter(dataset == "native") %>% sample_n(nrow), aes(x = PC1_PLM, y = PC2_PLM)) +
  geom_point(aes(x = PC1_PLM, y = PC2_PLM, colour = AbOrig_isotype), 
             size = 0.35, alpha = 1) +
  theme_bw()+
  xlab("PC1 (39%)") + ylab("PC2 (13%)") +
  scale_colour_manual(values = c("#7a0000", "#a80000", "#c10600", "#c45e05", "#c78219"),
                      name = "Isotype") +
  theme(axis.title = element_text(size = 12), legend.position = "top", legend.justification = "left",
        legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        text = element_text(family="Helvetica")) +
  guides(colour = guide_legend(override.aes = list(size=6)))

ggsave(p, file = "~/developability_project/plots/PCAs/integration_supporting_figure_PLM_by_isotype.png", height = 5.5, width = 5)


ggplot(DPL_PLM_PCAs_full_light %>% filter(dataset == "native"), aes(x = PC1_PLM, y = PC2_PLM))+
  stat_binhex(aes(x = PC1_PLM, y = PC2_PLM), bins = 100)+ 
  geom_point(data = DPL_PLM_PCAs_full_light %>% filter(dataset == "thera"),
             aes(x = PC1_PLM, y = PC2_PLM), 
             shape = 21,
             fill = "#5e3c99",
             colour = "white",
             size = 2, alpha = 0.7) +
  theme_bw()+
  theme(#legend.position = "top", legend.justification = "left",legend.key.width= unit(1.5, 'cm'), #run once for the legend
    legend.position = "none",
    text = element_text(family="Helvetica"),
    axis.title = element_text(size = 12)) +
  labs(fill = "Native antibody count") +
  xlab("PC1 (37%)") + ylab("PC2 (13%)")


ggsave("~/developability_project/plots/PCAs/integration_PCA_PLM_embeddings_native_and_mabs_hex_legend.pdf", height = 6, width = 6)
ggsave("~/developability_project/plots/PCAs/integration_PCA_PLM_embeddings_native_and_mabs_hex_light.png", height = 6, width = 9)
ggsave("~/developability_project/plots/PCAs/integration_PCA_PLM_embeddings_native_and_mabs_hex_light.pdf", height = 6, width = 9)




p = ggplot(DPL_PLM_PCAs_full_light %>% filter(dataset == "native"), aes(x = PC1_PLM, y = PC2_PLM))+
  stat_binhex(aes(x = PC1_PLM, y = PC2_PLM), bins = 100)+ 
  geom_point(data = DPL_PLM_PCAs_full_light %>% filter(dataset == "patented"),
             aes(x = PC1_PLM, y = PC2_PLM), colour = "#b2abd2", 
             size = 0.5, alpha = 0.5) +
  theme_bw()+
  theme(legend.position = "none",
        text = element_text(family="Helvetica"),
        axis.title = element_text(size = 12)) +
  labs(fill = "Native antibody count") +
  xlab("PC1 (37%)") + ylab("PC2 (13%)")


ggsave(p, file = "~/developability_project/plots/PCAs/integration_PCA_PLM_embeddings_native_and_pads_hex_light.png", height = 6, width = 9)
ggsave(p, file = "~/developability_project/plots/PCAs/integration_PCA_PLM_embeddings_native_and_pads_hex_light.pdf", height = 6, width = 9)

#check the role of isotype and IGHV in positioning in DPL space

nrow = DPL_PLM_PCAs_full_light %>% filter(dataset == "native") %>% nrow()
p = ggplot(DPL_PLM_PCAs_full_light %>% filter(dataset == "native") %>% 
             filter(AbOrig_isotype == "IgL")) +
  #sample_n(nrow), aes(x = PC1_DPL, y = PC2_DPL)) +
  geom_point(aes(x = PC1_DPL, y = PC2_DPL, colour = AbOrig_isotype), 
             size = 0.25, alpha = 1) +
  theme_bw()+
  xlab("PC1 (20%)") + ylab("PC2 (13%)") +
  #scale_colour_manual(values = c("#c9a847", "#cbcc82"),
  #name = "Isotype") +
  theme(axis.title = element_text(size = 12), legend.position = "top", legend.justification = "left",
        legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        text = element_text(family="Helvetica")) +
  guides(colour = guide_legend(override.aes = list(size=6)))

ggsave(p, file = "~/developability_project/plots/PCAs/integration_supporting_figure_by_isotype_light_DPL.png", height = 5.5, width = 5)
#ggsave(p, file = "~/developability_project/plots/PCAs/integration_supporting_figure_by_isotype.pdf", height = 5.5, width = 5)

DPL_PLM_PCAs_full_light %>% filter(dataset == "native") %>% select(AbOrig_v_gene_prefix) %>% drop_na() #852,326/854,408 --> 2K has no v-gene annotations?

p = ggplot(DPL_PLM_PCAs_full_light %>% filter(dataset == "native") %>% sample_n(nrow), aes(x = PC1_DPL, y = PC2_DPL)) +
  geom_point(aes(x = PC1_DPL, y = PC2_DPL, colour = AbOrig_v_gene_prefix), 
             size = 0.35, alpha = 1) +
  theme_bw()+
  xlab("PC1 (20%)") + ylab("PC2 (13%)") +
  scale_colour_brewer(palette = "RdYlBu", name = "IGHV family") +
  theme(axis.title = element_text(size = 12), legend.position = "top", legend.justification = "left",
        legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        text = element_text(family="Helvetica")) +
  guides(#colour = "none",
    colour = guide_legend(override.aes = list(size=6)))



ggsave(p, file = "~/developability_project/plots/PCAs/integration_supporting_figure_PLM_by_vgene_light.png", height = 5.5, width = 5)
# too many coplours for V genes, split the data to IgK and IgL 

p = ggplot(DPL_PLM_PCAs_full_light %>% filter(dataset == "native") %>%  
             filter(str_detect(AbOrig_v_gene_prefix, "IGL")) %>% #184,146 
             sample_n(184146 ), aes(x = PC1_DPL, y = PC2_DPL)) +
  geom_point(aes(x = PC1_DPL, y = PC2_DPL, colour = AbOrig_v_gene_prefix), 
             size = 0.35, alpha = 1) +
  theme_bw()+
  xlab("PC1 (20%)") + ylab("PC2 (13%)") +
  scale_colour_brewer(palette = "RdYlBu", name = "") +
  theme(axis.title = element_text(size = 12), legend.position = "top", legend.justification = "left",
        legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        text = element_text(family="Helvetica")) +
  guides(#colour = "none",
    colour = guide_legend(override.aes = list(size=6)))

ggsave(p, file = "~/developability_project/plots/PCAs/integration_supporting_figure_DPL_by_vgene_light_IgL.png", height = 5.5, width = 5)


p = ggplot(DPL_PLM_PCAs_full_light %>% filter(dataset == "native") %>%  
             filter(str_detect(AbOrig_v_gene_prefix, "IGK")) %>% #183,915 
             sample_n(183915), aes(x = PC1_DPL, y = PC2_DPL)) +
  geom_point(aes(x = PC1_DPL, y = PC2_DPL, colour = AbOrig_v_gene_prefix), 
             size = 0.35, alpha = 1) +
  theme_bw()+
  xlab("PC1 (20%)") + ylab("PC2 (13%)") +
  scale_colour_brewer(palette = "RdYlBu", name = "") +
  theme(axis.title = element_text(size = 12), legend.position = "top", legend.justification = "left",
        legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        text = element_text(family="Helvetica")) +
  guides(#colour = "none",
    colour = guide_legend(override.aes = list(size=6)))

ggsave(p, file = "~/developability_project/plots/PCAs/integration_supporting_figure_DPL_by_vgene_light_IgK.png", height = 5.5, width = 5)








#Figure 7C-------



#Matteo's 
read_non_native_predictability = function(file_input, dataset_name){
  output = read_csv(paste0("~/path/to/files/", file_input)) %>%
    select(starts_with(c("AbChain", "AbStruc"))) %>%
    rowid_to_column() %>%
    rename(repetition = rowid) %>%
    mutate(dataset = dataset_name) %>%
    select(repetition, dataset, everything())
  return(output)
}



PLM_predictability_humanised = read_non_native_predictability("best_regressors_humanised_scores.csv", "Kymouse") 
PLM_predictability_PAD = read_non_native_predictability("best_regressors_pad_scores.csv", "PAD") 
PLM_predictability_mAbs = read_non_native_predictability("best_regressors_thera_scores.csv", "mAbs")



#Jahn's


DPL_predictions_non_native_means = readRDS("~/developability_project/plotting_dataframes/from_jahn/predictability_july_2023/predictions_non_native.RDS") %>% 
  mutate(rsq = if_else(rsq < 0, 0, rsq)) %>% #make rsq values which are less than 0 == 0 
  group_by(dp_set_name, test_set, parameter) %>%
  summarise(mean_rsq = mean(rsq), sd_rsq = sd(rsq)) %>%
  mutate(level = case_when(str_detect(parameter, "AbChain") ~ "Sequence",
                           str_detect(parameter, "AbStruc") ~ "Structure"))




#data merging and curating

DPL_predictions_non_native_means_plotting = DPL_predictions_non_native_means %>% 
  filter(dp_set_name == "comb_parameters_mwds") %>%
  filter(test_set != "native") %>% 
  mutate(test_set = str_replace_all(test_set, "humanized", "Kymouse")) %>%
  mutate(test_set = str_replace_all(test_set, "patented", "PAD")) %>%
  mutate(test_set = str_replace_all(test_set, "therapeutic", "mAbs")) 

PLM_predictability_humanised = PLM_predictability_humanised %>%
  pivot_longer(-c(repetition,dataset), names_to = "parameter", values_to = "rsq") %>%
  filter(parameter %in% DPL_predictions_non_native_means_plotting$parameter) %>% 
  group_by(parameter,dataset) %>%
  summarise(mean_rsq = mean(rsq), sd_rsq = sd(rsq))

PLM_predictability_PAD  = PLM_predictability_PAD %>%
  pivot_longer(-c(repetition,dataset), names_to = "parameter", values_to = "rsq") %>%
  filter(parameter %in% DPL_predictions_non_native_means_plotting$parameter) %>% 
  group_by(parameter,dataset) %>%
  summarise(mean_rsq = mean(rsq), sd_rsq = sd(rsq))

PLM_predictability_mAbs = PLM_predictability_mAbs %>%
  pivot_longer(-c(repetition,dataset), names_to = "parameter", values_to = "rsq") %>%
  filter(parameter %in% DPL_predictions_non_native_means_plotting$parameter) %>% 
  group_by(parameter,dataset) %>%
  summarise(mean_rsq = mean(rsq), sd_rsq = sd(rsq))


PLM_predictions_non_native_means_plotting = rbind(PLM_predictability_humanised, PLM_predictability_PAD,PLM_predictability_mAbs) %>%
  mutate(level = case_when(str_detect(parameter, "AbChain") ~ "Sequence",
                           str_detect(parameter, "AbStruc") ~ "Structure")) %>%
  mutate(embedding = "PLM")


DPL_predictions_non_native_means_plotting = DPL_predictions_non_native_means_plotting %>%
  ungroup() %>%
  rename(dataset = test_set) %>%
  select(-dp_set_name) %>%
  mutate(embedding = "DPLs")


non_native_means_plotting = rbind(DPL_predictions_non_native_means_plotting,
                                  PLM_predictions_non_native_means_plotting)
save(non_native_means_plotting, file = "~/developability_project/environments_DataObjects/non_native_means_plotting.RData")




non_native_means_plotting = non_native_means_plotting %>%
  rename(Embedding = embedding)

non_native_means_plotting$dataset = factor(non_native_means_plotting$dataset, levels = c("Kymouse", "PAD", "mAbs"))

non_native_means_plotting_text = non_native_means_plotting %>%
  group_by(dataset, Embedding, level) %>%
  summarise(MD = median(mean_rsq),
            mean = mean(mean_rsq))

ggplot(non_native_means_plotting, aes(x = dataset, y = mean_rsq, fill = Embedding)) +
  facet_grid(rows = vars(level)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  geom_jitter(alpha=0.6, pch = 21, size = 1.5, stroke = 0, aes(fill = Embedding), position = position_jitterdodge()) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02")) +
  theme_bw() +
  geom_text(data = non_native_means_plotting_text,
            aes(x = dataset, y = -0.1, label = round(mean, digits = 2), color = Embedding), position = position_dodge(width = 0.9),
            show.legend = F)+
  theme(strip.background = element_rect(fill = "white"), strip.text = element_text(size = 12),
        legend.position = "top", legend.justification = "left",
        axis.text = element_text(size = 10),
        text = element_text(family="Helvetica")) +
  scale_colour_manual(values = c("#1b9e77", "#d95f02")) +
  ylab(expression(Mean~R^2)) + xlab("")

ggsave("~/developability_project/plots/predictability/non_native_from_native_training_DPLs_PLM_MWDSonly.png", width = 5, height = 8)
ggsave("~/developability_project/plots/predictability/non_native_from_native_training_DPLs_PLM_MWDSonly.pdf", width = 5, height = 8)





