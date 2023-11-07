#Author: Habib Bashour
library(ComplexHeatmap)
library(reshape2)
library(rstatix)
library(ggpubr)
library(cowplot)

#correlation matrix production------
#sequence-only DPs
AbChain_native_developability_final = read_csv("~/data/native/AbChain_native_developability_final.csv")

isotypes = c("IgG", "IgA", "IgM", "IgD" ,"IgE" ,"IgK" ,"IgL" ,"IgG", "IgM" ,"IgK", "IgL")
species  = c("human" ,"human", "human", "human" ,"human", "human", "human" ,"mouse", "mouse" ,"mouse","mouse")


for (i in 1:length(isotypes)) {
  x =  isotypes[i]
  y = species[i]
  data = AbChain_native_developability_final %>%
    filter(isotype == x & species == y) %>%
    select(-c(rowid, isotype ,species, chain, starts_with("aaSeq"), database, v_gene, j_gene)) 
  AbChain_parameter_cormatrix_fv =  round(cor(data, method = "pearson", use = "pairwise.complete.obs"),2)
  print(i)
  save(AbChain_parameter_cormatrix_fv, file = paste0("~/developability_project/environments_DataObjects/parameter_correlation_matrices_pearson/fv/AbChain_parameter_cormatrix_fv_",isotypes[i],"_",species[i],".RData"))

}

#structure-only DPs
AbStruc_native_developability_final = read_csv("~/data/native/AbStruc_native_developability_final.csv")

for (i in 1:length(isotypes)) {
  x =  isotypes[i]
  y = species[i]
  data = AbStruc_native_developability_final %>%
    filter(isotype == x & species == y) %>%
    select(-c(rowid, isotype ,species, chain, starts_with("aaSeq"), database, v_gene, j_gene)) 
  
  AbStruc_parameter_cormatrix_fv =  round(cor(data, method = "pearson", use = "pairwise.complete.obs"),2)

  print(i)
  save(AbStruc_parameter_cormatrix_fv, file = paste0("~/developability_project/environments_DataObjects/parameter_correlation_matrices_pearson/fv/AbStruc_parameter_cormatrix_fv_",isotypes[i],"_",species[i],".RData"))
}

#both levels (sequence and structure DPs)

fv_developability =  AbChain_native_developability_final%>% full_join(AbStruc_native_developability_final)

for (i in 1:length(isotypes)) {
  x =  isotypes[i]
  y = species[i]
  data = fv_developability %>%
    filter(isotype == x & species == y) %>%
    select(-c(rowid, isotype ,species, chain, starts_with("aaSeq"), database, v_gene, j_gene)) 

    AbChStruc_parameter_cormatrix_fv =  round(cor(data, method = "pearson", use = "pairwise.complete.obs"),2)
    print(i)
  
  save(AbChStruc_parameter_cormatrix_fv, file = paste0("~/developability_project/environments_DataObjects/parameter_correlation_matrices_pearson/fv/AbChStruc_parameter_cormatrix_fv_",isotypes[i],"_",species[i],".RData"))
}





#produce the overall pairwise DP correlation stats for Figure 2A----------

setwd("~/developability_project/environments_DataObjects/parameter_correlation_matrices_pearson/fv/") #the directory where all the correlation matrices are saved
list_RData = dir(pattern = "*.RData")

correlations_stats = tibble(kind = character(), category = character(),
                            perason_cor = double(), Var1 = factor(),Var2 = factor())

for (d in 1:length(list_RData)){
  matrix = loadRData(paste0("~/developability_project/environments_DataObjects/parameter_correlation_matrices_pearson/fv/", list_RData[d]))
  diag(matrix) <- NA
  matrix[lower.tri(matrix)] <- NA
  cat = str_remove_all(list_RData[d], "Ab.*_parameter_cormatrix_|.RData")
  kinD = str_remove_all(list_RData[d], "_param.*")
  single_correlations_stats = matrix %>%
    melt() %>%
    as_tibble() %>%
    mutate(value = abs(value)) %>%
    rename(perason_cor = value) %>%
    #ggplot(aes(x = perason_cor)) +
    #geom_density()
    mutate(category = cat, kind = kinD) %>%
    mutate(kind = str_replace_all(kind, "AbChain", "sequence"),
           kind = str_replace_all(kind, "AbStruc", "structure"),
           kind = str_replace_all(kind, "AbChStruc", "sequence_and_structure"))
  correlations_stats = rbind(correlations_stats,single_correlations_stats )
  
}



correlations_stats = correlations_stats %>%
  filter(kind != "sequence_and_structure") %>%
  mutate(category = str_remove_all(category, "fv_")) %>%
  separate(category, c("isotype", "species"), sep = "_")
mutate(species = str_to_title(species))

isotype_plotting_order = c("IgD","IgM","IgG","IgA","IgE","IgK","IgL")


correlations_stats$isotype = factor(correlations_stats$isotype, levels =isotype_plotting_order)
correlations_stats = correlations_stats %>%
  mutate(species = str_to_title(species))


correlations_stats_text = correlations_stats %>%
  group_by(isotype, species, kind) %>%
  summarize(MD = median(perason_cor))

correlations_stats_text$isotype = factor(correlations_stats_text$isotype, levels =isotype_plotting_order)

#perform statistical test between sequence and structure DP pairwise correlation

correlations_stats_test = correlations_stats %>%
  group_by(species, isotype) %>%
  wilcox_test(perason_cor ~ kind, comparisons = list(c("sequence", "structure"))) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() 


save(correlations_stats, file = "~/developability_project/environments_DataObjects/correlations_stats_20221031.RData")
load("~/developability_project/environments_DataObjects/correlations_stats_20221031.RData")

ggplot(correlations_stats, aes(x = isotype, y = perason_cor)) +
  facet_grid(cols = vars(species), scales = "free_x", space = "free_x") +
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"), strip.text = element_text(size = 10),
        legend.position = "top", legend.justification="left") +
  scale_fill_manual(values = c("#9400D4", "#98a8ea")) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.9), width = 0.7, aes(fill = kind), alpha = 0.3) +
  geom_point(aes(fill = kind, colour = kind), position = position_jitterdodge(jitter.height = 0.01, jitter.width = 0.3, dodge.width = 0.9),
             pch = 21, size = 0.7, stroke = 0.1, alpha = 0.7) +
  scale_color_manual(values = c("#9400D4", "#98a8ea")) +
  labs(fill = "Level") +
  geom_text(data = correlations_stats_text, aes(x = isotype, y = MD+0.02, label = MD, group = kind), colour = "white", position = position_dodge(0.9), size = 3) +
  ylab("Absolute Pearson correlation coefficient") + xlab("Isotype") 




ggsave2("~/developability_project/plots/parameter_correlations/updated_paramters/fv/absolute_paramter_correlation_boxplot_distribution.png", width = 9, height = 5.5)
ggsave2("~/developability_project/plots/parameter_correlations/updated_paramters/fv/absolute_paramter_correlation_boxplot_distribution.pdf", width = 9, height = 5.5)



#ABC-EDA algorithm output processing, provided by Jahn Zhong October 2022 ------

ABC_EDA_output_fv_updated = read.csv("~/developability_project/environments_DataObjects/dominating_sets/ABC_EDA/updated_paramters/sequence_and_structure_ds.csv", header = TRUE) %>%
  as_tibble() %>%
  separate(subset_threshold, c("kind", "region","isotype", "species", "thresh"), sep = "_") %>%
  mutate(thresh = as.double(thresh)) %>%
  mutate(cat = paste0(isotype, "_", species)) %>%
  mutate(dominating_set = str_remove_all(dominating_set,"'")) %>%
  mutate(dominating_set = str_remove_all(dominating_set," ")) %>%
  mutate(isolated_nodes = str_remove_all(isolated_nodes,"'")) %>%
  mutate(isolated_nodes = str_remove_all(isolated_nodes," ")) %>%
  dplyr::rename(ABC_EDA = dominating_set) %>%
  mutate(doublets = str_remove_all(doublets, "'")) %>%
  mutate(doublets = str_remove_all(doublets, " ")) %>%
  mutate(doublets = str_remove_all(doublets,"\\{")) %>%
  mutate(doublets = str_remove_all(doublets,"\\}"))


ABC_EDA_output_fv_updated = ABC_EDA_output_fv_updated %>%
  mutate(cat = paste0(isotype, "_", species))





#heatmap annotation function-------

get_parameter_correlation_annotation_df = function(input_matrix){
  as_tibble(colnames(input_matrix)) %>%
    dplyr::rename(measure = value) %>%
    filter(!measure %in% c("rowid", "cloneId", "isotype", "database", "species","chain","v_gene", "j_gene")) %>%
    mutate(measure = str_replace_all(measure, "lenght", "length")) %>%
    mutate(group = case_when(str_detect(measure, "_mw") | str_detect( measure, "_length") | str_detect( measure, "_av_resdiue_weight") ~ "molecular",
                             str_detect(measure, "content") | str_detect( measure, "free_cys") | str_detect(measure, "cys_bridges") ~ "aa composition",
                             str_detect(measure, "extcoef") ~ "photochemical",
                             str_detect(measure, "index") ~ "stability",
                             str_detect(measure, "charge") | str_detect( measure, "pI") | str_detect( measure, "hydrophobicity") | str_detect( measure, "hmom") | str_detect( measure, "solubility") ~ "electrochemical",
                             str_detect(measure, "rank") | str_detect( measure, "binders") | str_detect( measure, "span") ~ "immunogenicity",
                             str_detect(measure, "angle") | str_detect( measure, "beta") | str_detect( measure, "helices") | str_detect( measure, "loops") ~ "secondary structure",
                             str_detect(measure, "energy") ~ "thermodynamic",
                             str_detect(measure, "lauer") ~ "druggability",
                             str_detect(measure, "sasa") ~ "solvent accessibility",
                             TRUE ~ "interactions")) %>%
    mutate(type = case_when(str_detect(measure, "AbChain") ~ "sequence",
                            str_detect(measure, "AbStruc") ~ "structure")) %>%
    mutate(region = case_when(str_detect(measure, "fr1") ~ "FR1",str_detect(measure, "fr2") ~ "FR2",str_detect(measure, "fr3") ~ "FR3",str_detect(measure, "fr4") ~ "FR4",
                              str_detect(measure, "cdr1") ~ "CDR1",str_detect(measure, "cdr2") ~ "CDR2",str_detect(measure, "cdr3") ~ "CDR3",
                              str_detect(measure, "AbChain") | str_detect( measure, "AbStruc") ~ "Fv"))
}

#heatmap production for all isotype/species combos with annotation of ABC-EDA algorithm output for all PEarson correlation thresholds (0.1,0.9.0.1)----


setwd("~/developability_project/environments_DataObjects/parameter_correlation_matrices_pearson/fv/") #the directory where all the correlation matrices are saved
list_RData = dir(pattern = "*.RData")
#plot the ds_anotated heatmaps
for (j in 1:length(list_RData)){
  
  file_name = str_replace(list_RData[j], "cormatrix", "corheatmap")
  file_name = str_remove(file_name, ".RData")
  plot_name = str_remove(file_name, "Ab.+_parameter_corheatmap_")
  plot_name = str_remove(plot_name, ".pdf")
  category = str_remove(plot_name, "fv_")
  kinD = str_remove(file_name, "_parameter.*")
  matrix = loadRData(paste0("~/developability_project/environments_DataObjects/parameter_correlation_matrices_pearson/fv/", list_RData[j]))
  diag(matrix) = NA
  annotation_df = get_parameter_correlation_annotation_df(matrix)
  for(i in 1:length(thresholds)){
    ABC_EDA = ABC_EDA_output_fv_updated %>% filter(cat == category & thresh == thresholds[i] & kind == kinD) %>% pull(ABC_EDA) %>% unlist() %>% str_split(",") %>% unlist()
    isolated_nodes = ABC_EDA_output_fv_updated %>% filter(cat == category & thresh == thresholds[i] & kind == kinD) %>% pull(isolated_nodes) %>% unlist() %>% str_split(",") %>% unlist()
    doublets = ABC_EDA_output_fv_updated %>% filter(cat == category & thresh == thresholds[i]  & kind == kinD) %>% pull(doublets) %>% unlist() %>% str_split(",") %>% unlist()
    annotation_df = annotation_df %>%
      mutate(dominance = case_when(#measure %in% ABC_EDA ~ "MWDS",
        measure %in% ABC_EDA ~ "dominant",
        measure %in% doublets ~ "doublet element",
        measure %in% isolated_nodes ~ "isolated node",
        TRUE ~ "redundant"))
    #annotation objects-----
    group_annotation_bottom = HeatmapAnnotation(`Physicochemical property` = annotation_df$group,
                                                col = list(`Physicochemical property` = c("aa composition" = "#543005", "druggability" = "#8c510a",
                                                                                          "electrochemical" = "#bf812d", "interactions" = "#dfc27d", "molecular" = "#f6e8c3",
                                                                                          "immunogenicity" = "#f5f5f5","photochemical" = "#c7eae5", "secondary structure" = "#80cdc1",
                                                                                          "solvent accessibility" = "#35978f", "stability" = "#01665e","thermodynamic" = "chocolate")))#,show_legend = FALSE, show_annotation_name = FALSE)
    
    type_annotation_bottom = HeatmapAnnotation(Level = annotation_df$type, col = list(Level = c("sequence" = "#9400D4", "structure" = "#98a8ea" )))#,show_legend = FALSE, show_annotation_name = FALSE)
    
    
    region_annotation_bottom =  HeatmapAnnotation(region = annotation_df$region ,col = list(region = c("Fv" = "salmon", "CDR1" = "orange", "CDR2" = "orange3", "CDR3" = "orange4",
                                                                                                       "FR1" = "limegreen", "FR2" = "seagreen1","FR3" = "palegreen2","FR4" = "seagreen4")),show_legend = FALSE, show_annotation_name = FALSE)
    
    dominance_annotation_bottom =  HeatmapAnnotation(Dominance = annotation_df$dominance, col = list(Dominance = c(#"MWDS" = "forestgreen",
      "dominant" = "#b2182b",
      "isolated node" = "#67a9cf","doublet element" = "#d1e5f0",
      "redundant" = "#f4a582")))#,
    #annotation_legend_param = list(labels = c("dominant", "isolated node", "doublet element", "redundant")))#,show_legend = FALSE, show_annotation_name = FALSE)
    
    
    group_annotation_right = rowAnnotation(`Physicochemical property` = annotation_df$group,
                                           col = list(`Physicochemical property` = c("aa composition" = "#543005", "druggability" = "#8c510a",
                                                                                     "electrochemical" = "#bf812d", "interactions" = "#dfc27d", "molecular" = "#f6e8c3",
                                                                                     "immunogenicity" = "#f5f5f5","photochemical" = "#c7eae5", "secondary structure" = "#80cdc1",
                                                                                     "solvent accessibility" = "#35978f", "stability" = "#01665e","thermodynamic" = "chocolate")),show_legend = FALSE, show_annotation_name = FALSE)
    
    type_annotation_right = rowAnnotation(Level = annotation_df$type, col = list(Level = c("sequence" = "#9400D4", "structure" = "#98a8ea")),show_legend = FALSE, show_annotation_name = FALSE)
    
    
    region_annotation_right =  rowAnnotation(region = annotation_df$region ,col = list(region = c("Fv" = "salmon", "CDR1" = "orange", "CDR2" = "orange3", "CDR3" = "orange4",
                                                                                                  "FR1" = "limegreen", "FR2" = "seagreen1","FR3" = "palegreen2","FR4" = "seagreen4")))
    dominance_annotation_right =  rowAnnotation(Dominance = annotation_df$dominance, col = list(Dominance = c("MWDS" = "forestgreen",
                                                                                                              "dominant" = "#b2182b",
                                                                                                              "isolated node" = "#67a9cf",
                                                                                                              "doublet element" = "#d1e5f0",
                                                                                                              "redundant" = "#f4a582")),show_legend = FALSE, show_annotation_name = FALSE)
    
    
    pdf(file= paste0("~/developability_project/plots/parameter_correlations/updated_paramters/fv/ds_annotated/", file_name, "_", thresholds[i],".pdf"), width=15, height=11)
    #plotting----
    hmap = Heatmap(matrix, column_title = "" , name = "Pearson correlation",
                   show_column_dend = F, show_row_dend = F,
                   #row_order = order, column_order = order,
                   #rect_gp = gpar(col = "black", lwd = 1),
                   border = F, na_col = "white",
                   column_title_gp = gpar(fill = "white", col = "black", border = "white"),
                   top_annotation = c(type_annotation_bottom,group_annotation_bottom,dominance_annotation_bottom),
                   right_annotation = c(dominance_annotation_right,group_annotation_right, type_annotation_right),
                   #top_annotation = c(region_annotation_bottom),
                   
                   row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize = 5),
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf("%.1f", matrix[i, j]), x, y, gp = gpar(fontsize =2.5, col = "black"))
                   })
    
    draw(hmap, merge_legend = TRUE)
    dev.off()
    
  }
} 
#total of 297 heatmaps (11X9X3), among those we used the IgG human subset on both levels for 0.6 threshold (AbChStruc_parameter_corheatmap_fv_IgE_human_0.6.pdf) for Figure2B
#the remaining heatmaps (also both levels) for the rest of isotypes were used for Supp. Figures 2 and 3







#Supp Figure 4A (the count of isolated noeds and dominant parameters)------ 

ABC_EDA_output_plotting = ABC_EDA_output_fv_updated %>%
  filter(kind =="AbChStruc") %>%
  mutate(#MWDS_count = str_count(ABC_EDA,",") +1,
    #doublet_count = str_count(doublets,",") +1,
    #isolated_nodes_count = str_count(isolated_nodes,",") +1,
    MWDS_AbChain = str_count(ABC_EDA,"AbChain"),
    MWDS_AbStruc = str_count(ABC_EDA,"AbStruc"),
    doublet_AbChain = str_count(doublets,"AbChain"),
    doublet_AbStruc = str_count(doublets,"AbStruc"),
    isolated_nodes_AbChain = str_count(isolated_nodes,"AbChain"),
    isolated_nodes_AbStruc = str_count(isolated_nodes,"AbStruc"),
    MWDS_count = MWDS_AbChain + MWDS_AbStruc ,
    doublet_count = doublet_AbChain + doublet_AbStruc,
    isolated_nodes_count = isolated_nodes_AbChain + isolated_nodes_AbStruc) %>%
  select(-c(ABC_EDA ,doublets, isolated_nodes)) %>%
  select(-c(isotype, species)) %>%
  pivot_longer(-c(kind,region,thresh,cat), names_to ="identifier", values_to = "paramter_count") %>%
  mutate(identifier_type = case_when(str_detect(identifier, "MWDS") ~ "dominant",
                                     str_detect(identifier, "doublet") ~ "doublet",
                                     str_detect(identifier, "isolated_node") ~ "isolated node")) %>%
  mutate(level = case_when(str_detect(identifier, "AbChain") ~ "sequence",
                           str_detect(identifier, "AbStruc") ~ "structure"))
save(ABC_EDA_output_plotting, file = "~/developability_project/environments_DataObjects/ABC_EDA_output_plotting.RData")

ABC_EDA_output_plotting$cat = factor(ABC_EDA_output_plotting$cat, levels = c("IgD_human", "IgM_human", "IgG_human", "IgA_human", "IgE_human",
                                                                             "IgK_human", "IgL_human", "IgM_mouse", "IgG_mouse", "IgK_mouse", "IgL_mouse"))
ABC_EDA_output_plotting = ABC_EDA_output_plotting %>%
  rename(Level = level)

save(ABC_EDA_output_plotting, file = "~/developability_project/environments_DataObjects/ABC_EDA_output_plotting.RData")

ABC_EDA_performance_evaluation_isolaed_nodes = ABC_EDA_output_plotting %>%
  filter(identifier_type == "isolated node") %>%
  mutate(paramter_proportion = case_when(Level == "sequence" ~ (paramter_count/40)*100, 
                                         Level == "structure" ~ (paramter_count/46)*100)) 

ggplot(ABC_EDA_performance_evaluation_isolaed_nodes, aes(x = thresh, y = paramter_proportion)) +
  geom_col(position = position_dodge(width = 0.1), aes(fill = Level)) +
  facet_wrap(vars(cat)) +
  theme_bw() +
  geom_text(data = ABC_EDA_performance_evaluation_isolaed_nodes, aes(x = thresh, y = paramter_proportion+1.5, label = round(paramter_proportion,digits =1), group = Level), position = position_dodge(width = 0.1), size = 2.5) +
  scale_x_continuous(breaks = seq(0.1,0.9,0.1)) +
  xlab("Pearson correlation coeffecient threshold") +
  ylab("Proportion to the total number of DPs within the same level (%)") +
  scale_fill_manual(values= c("#9400D4","#98a8ea")) +
  #ggtitle("Dominance class: isolated nodes", subtitle = "starting number of DPs = 88") +
  theme(strip.background = element_rect(fill = "white"), strip.text = element_text(size = 9), legend.position = "top",
        legend.justification = "left",
        axis.title = element_text(size = 14), legend.text = element_text(size = 12))

ggsave("~/developability_project/plots/ABC_EDA_paramter_optimisation/ABC_EDA_performance_evaluation_isolaed_nodes_proportion.png" , width = 15, height = 9.5)
ggsave("~/developability_project/plots/ABC_EDA_paramter_optimisation/ABC_EDA_performance_evaluation_isolaed_nodes_proportion.pdf" , width = 15, height = 9.5)
save(ABC_EDA_performance_evaluation_isolaed_nodes, file = "~/developability_project/environments_DataObjects/ABC_EDA_performance_evaluation_isolaed_nodes.RData")

#for the text:



ABC_EDA_performance_evaluation_dominants = ABC_EDA_output_plotting %>%
  drop_na() %>%
  filter(identifier_type == "dominant") %>%
  mutate(paramter_proportion = case_when(Level == "sequence" ~ (paramter_count/40)*100,
                                         Level == "structure" ~ (paramter_count/46)*100))

ggplot(ABC_EDA_performance_evaluation_dominants, aes(x = thresh, y = paramter_count)) +
  geom_col(position = position_dodge(width = 0.1), aes(fill = Level)) +
  facet_wrap(vars(cat)) +
  theme_bw() +
  geom_text(data = ABC_EDA_performance_evaluation_dominants, aes(x = thresh, y = paramter_count+1, label = paramter_count, group = Level), position = position_dodge(width = 0.1), size = 2.5) +
  scale_x_continuous(breaks = seq(0.1,0.9,0.1)) +
  xlab("Pearson correlation coeffecient threshold") +
  ylab("DP count") +
  scale_fill_manual(values= c("#9400D4","#98a8ea")) +
  #ggtitle("Dominance class: dominant parameters", subtitle = "starting number of DPs = 88") +
  theme(strip.background = element_rect(fill = "white"), strip.text = element_text(size = 9), legend.position = "top",
        legend.justification = "left",
        axis.title = element_text(size = 14), legend.text = element_text(size = 12))

ggsave("~/developability_project/plots/ABC_EDA_paramter_optimisation/ABC_EDA_performance_evaluation_dominants_count.png" , width = 15, height = 9.5)
ggsave("~/developability_project/plots/ABC_EDA_paramter_optimisation/ABC_EDA_performance_evaluation_dominants_count.pdf" , width = 15, height = 9.5)



ggplot(ABC_EDA_performance_evaluation_dominants, aes(x = thresh, y = paramter_proportion)) +
  geom_col(position = position_dodge(width = 0.1), aes(fill = Level)) +
  facet_wrap(vars(cat)) +
  theme_bw() +
  geom_text(data =ABC_EDA_performance_evaluation_dominants, aes(x = thresh, y = paramter_proportion+1.5, label = round(paramter_proportion,digits =1), group = Level), position = position_dodge(width = 0.1), size = 2.5) +
  scale_x_continuous(breaks = seq(0.1,0.9,0.1)) +
  xlab("Pearson correlation coeffecient threshold") +
  ylab("Proportion to the total number of DPs within the same level (%)") +
  scale_fill_manual(values= c("#9400D4","#98a8ea")) +
  #ggtitle("Dominance class: dominant parameters", subtitle = "starting number of DPs = 88") +
  theme(strip.background = element_rect(fill = "white"), strip.text = element_text(size = 9), legend.position = "top",
        legend.justification = "left",
        axis.title = element_text(size = 14), legend.text = element_text(size = 12))




  








