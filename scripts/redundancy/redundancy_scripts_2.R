#Author: Habib Bashour

library(ComplexHeatmap)
library(factoextra)
library(dendextend)
library(reticulate)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
} # a function to load and rename .RData objects 

#This script to produce the figures shown in main Figure 3 and its supporting Figures 
# object ABC_EDA_output_fv_updated is provided as sequence_and_structure_ds.csv and processed in redundancy_scripts.R
#Panel A - sequence DPs-----

ABC_EDA_output_MWDS = ABC_EDA_output_fv_updated %>%
  filter(thresh == 0.6) %>%
  mutate(MWDS = paste0(ABC_EDA, ",", doublets, ",", isolated_nodes))


sequence_ABC_EDA_output_MWDS = ABC_EDA_output_MWDS %>%
  filter(thresh == 0.6 & kind == "AbChain")


for (i in 1:nrow(sequence_ABC_EDA_output_MWDS)){
  category = pluck(sequence_ABC_EDA_output_MWDS, "cat", i)
  sequence_MWDS_0.6_list[[category]] = pluck(sequence_ABC_EDA_output_MWDS, "MWDS", i) %>% str_split(",") %>% unlist() %>% c()
}




get_upset_plot = function(input_comb_matrix, ss = set_size(input_comb_matrix),cs = comb_size(input_comb_matrix), order_input){
  ht = UpSet(input_comb_matrix,
             #set_order = order(ss),
             set_order = order_input,
             comb_order = order(comb_degree(input_comb_matrix), -cs),
             top_annotation = HeatmapAnnotation(
               "MWDS intersection" = anno_barplot(cs,
                                                  ylim = c(0, max(cs)*1.1),
                                                  border = FALSE,
                                                  gp = gpar(fill = "black"),
                                                  height = unit(4, "cm")),
               annotation_name_side = "left",
               annotation_name_rot = 90),
             left_annotation = rowAnnotation(
               #"MWDS paramter count" = anno_barplot(-ss),
               set_name = anno_text(set_name(input_comb_matrix),
                                    location = 0.5,
                                    just = "center",
                                    width = max_text_width(set_name(input_comb_matrix)) + unit(4, "mm"))
             ),
             right_annotation = NULL,
             show_row_names = FALSE)
  ht = draw(ht)
  od = column_order(ht)
  decorate_annotation("MWDS intersection", {
    grid.text(cs[od], x = seq_along(cs)-0.17, y = unit(cs[od], "native") + unit(2, "pt"),
              default.units = "native", just = c("left", "bottom"),
              gp = gpar(fontsize = 11.25, col = "#404040"))
    #, rot = 45)
  })
}


heavy_human_comb_mat = make_comb_mat(IgA_human = sequence_MWDS_0.6_list$IgA_human,
                                     IgD_human = sequence_MWDS_0.6_list$IgD_human,
                                     IgE_human = sequence_MWDS_0.6_list$IgE_human,
                                     IgG_human = sequence_MWDS_0.6_list$IgG_human,
                                     IgM_human = sequence_MWDS_0.6_list$IgM_human,
                                     mode = "intersect")


light_comb_mat = make_comb_mat(IgK_mouse = sequence_MWDS_0.6_list$IgK_mouse,
                               IgL_mouse = sequence_MWDS_0.6_list$IgL_mouse,
                               IgK_human = sequence_MWDS_0.6_list$IgK_human,
                               IgL_human = sequence_MWDS_0.6_list$IgL_human,
                               mode = "intersect")

heavy_human_mouse_comb_mat = make_comb_mat(IgG_mouse = sequence_MWDS_0.6_list$IgG_mouse,
                                           IgM_mouse = sequence_MWDS_0.6_list$IgM_mouse,
                                           IgG_human = sequence_MWDS_0.6_list$IgG_human,
                                           IgM_human = sequence_MWDS_0.6_list$IgM_human,
                                           mode = "intersect")


heavy_human_order = c("IgD_human", "IgM_human", "IgG_human", "IgA_human", "IgE_human")
pdf("~/developability_project/plots/parameter_overlap/upset_plots/heavy_human_seqparamters_upsetplot.pdf", width = 10, height = 3)
get_upset_plot(heavy_human_comb_mat, order_input = heavy_human_order) #Supp. Figure 4B top
dev.off()

heavy_human_mouse_order = c( "IgM_human", "IgG_human", "IgM_mouse", "IgG_mouse")
pdf("~/developability_project/plots/parameter_overlap/upset_plots/heavy_human_mouse_seqparamters_upsetplot.pdf", width = 8, height = 3)
get_upset_plot(heavy_human_mouse_comb_mat, order_input = heavy_human_mouse_order) # Figure 3, Panel A top left
dev.off()

light_human_mouse_order = c( "IgK_human", "IgL_human", "IgK_mouse", "IgL_mouse")
pdf("~/developability_project/plots/parameter_overlap/upset_plots/light_human_mouse_seqparamters_upsetplot.pdf", width = 8, height = 3)
get_upset_plot(light_comb_mat, order_input = light_human_mouse_order) #Figure 3, Panel A bottom left
dev.off()


#Panel A - structure DPs-----

structure_ABC_EDA_output_MWDS = ABC_EDA_output_MWDS %>%
  filter(thresh == 0.6 & kind == "AbStruc")


for (i in 1:nrow(structure_ABC_EDA_output_MWDS)){
  category = pluck(structure_ABC_EDA_output_MWDS, "cat", i)
  structure_MWDS_0.6_list[[category]] = pluck(structure_ABC_EDA_output_MWDS, "MWDS", i) %>% str_split(",") %>% unlist() %>% c()
}


heavy_human_mouse_comb_mat_struc = make_comb_mat(IgG_mouse = structure_MWDS_0.6_list$IgG_mouse,
                                                 IgM_mouse = structure_MWDS_0.6_list$IgM_mouse,
                                                 IgG_human = structure_MWDS_0.6_list$IgG_human,
                                                 IgM_human = structure_MWDS_0.6_list$IgM_human,
                                                 mode = "intersect")
light_comb_mat_struc =  make_comb_mat(IgK_mouse = structure_MWDS_0.6_list$IgK_mouse,
                                      IgL_mouse = structure_MWDS_0.6_list$IgL_mouse,
                                      IgK_human = structure_MWDS_0.6_list$IgK_human,
                                      IgL_human = structure_MWDS_0.6_list$IgL_human,
                                      mode = "intersect")

heavy_human_struc_comb_mat = make_comb_mat(IgA_human = structure_MWDS_0.6_list$IgA_human,
                                           IgD_human = structure_MWDS_0.6_list$IgD_human,
                                           IgE_human = structure_MWDS_0.6_list$IgE_human,
                                           IgG_human = structure_MWDS_0.6_list$IgG_human,
                                           IgM_human = structure_MWDS_0.6_list$IgM_human,
                                           mode = "intersect")



pdf("~/developability_project/plots/parameter_overlap/upset_plots/heavy_human_mouse_strucparamters_upsetplot.pdf", width = 8, height = 3)
get_upset_plot(heavy_human_mouse_comb_mat_struc, order_input = heavy_human_mouse_order) # Figure 3, Panel A top right
dev.off()
pdf("~/developability_project/plots/parameter_overlap/upset_plots/light_human_mouse_strucparamters_upsetplot.pdf", width = 8, height = 3)
get_upset_plot(light_comb_mat_struc, order_input = light_human_mouse_order) # Figure 3, Panel A bottom right
dev.off()

pdf("~/developability_project/plots/parameter_overlap/upset_plots/heavy_human_strucparamters_upsetplot.pdf", width = 10, height = 3)
get_upset_plot(heavy_human_struc_comb_mat, order_input = heavy_human_order) #Supp. Figure 4B bottom 
dev.off()






#Panel B----

scale_manual_fun_dist_optimised = function(x){
  (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
}

my_get_dist = function(x, method = "euclidean", stand = FALSE){
  if (stand) 
    x <- x %>% mutate_at(scale_manual_fun_dist_optimised, .vars = vars(-c(Var1,Var2))) #this is required 
  x <- x %>% select(-c(Var1,Var2))
  if (method %in% c("pearson", "spearman", "kendall")) {
    res.cor <- stats::cor(x, method = method, use = "pairwise.complete.obs") 
    res.dist <- stats::as.dist(1 - res.cor)
  }
  else res.dist <- stats::dist(x, method = method, ...)
  res.dist
  
}

setwd("~/developability_project/environments_DataObjects/parameter_correlation_matrices_pearson/fv/") #where the correlation matrices are saved as R objects

list_RData = dir(pattern = "AbChain") #run (1)
 

#run the lines 181 - 190 for i = 1 until drop_na() then run the following
ordered_cor_df = sample_pairwise_data %>% select(-value)

for(i in 1:length(list_RData)){
  matrix = loadRData(paste0("~/developability_project/environments_DataObjects/parameter_correlation_matrices_pearson/fv/",list_RData[i]))
  file_name = str_remove_all(list_RData[i], "AbChain_parameter_cormatrix_fv_|.RData")
  #file_name = str_remove_all(list_RData[i], "AbStruc_parameter_cormatrix_fv_|.RData") run (2)
  #file_name = str_remove_all(list_RData[i], "AbChStruc_parameter_cormatrix_fv_|.RData") run (3)
  diag(matrix) = NA
  matrix[lower.tri(matrix)] = NA
  sample_pairwise_data = matrix %>%
    melt() %>%
    as_tibble() %>%
    drop_na()
  sample_pairwise_data = setNames(sample_pairwise_data, c("Var1", "Var2",file_name))
  ordered_cor_df = full_join(ordered_cor_df,sample_pairwise_data) #full_join so we do not lose columns when a matrix is missing the pair
}

#rerun the loop for each if the following lines: 

ordered_cor_df_AbChain = ordered_cor_df #after run (1)
#rerun the loop for each if the following lines:
list_RData = dir(pattern = "AbStruc") #run (2)
#ordered_cor_df_AbStruc = ordered_cor_df #after run (2)


AbChain_cor_dist = my_get_dist(ordered_cor_df_AbChain, method = "pearson", stand = TRUE)
AbStruc_cor_dist = my_get_dist(ordered_cor_df_AbStruc, method = "pearson", stand = TRUE)

AbChain_dend <- AbChain_cor_dist %>% hclust %>% as.dendrogram
AbStruc_dend <- AbStruc_cor_dist %>% hclust %>% as.dendrogram

#Figure 3B, left-----
pdf("~/developability_project/plots/cor_correlation/dendrograms/pearson_distance_pairwise_native_AbChain_uncoloured.pdf", width = 7, height = 5.5)
AbChain_dend %>%
  set("branches_lwd", 2) %>%
  plot(#main = "cor-cor sequence parameters",
    ylab = "correlation distance")
dev.off()


#Figure 3B, right----
pdf("~/developability_project/plots/cor_correlation/dendrograms/pearson_distance_pairwise_native_AbStruc_uncoloured.pdf", width = 7, height = 5.5)
AbChStruc_dend %>%
  set("branches_lwd", 2) %>%
  plot(#main = "cor-cor structure parameters",
    ylab = "correlation distance")
dev.off()




#Panel C----




query_4_data = read_csv("~/developability_project/plotting_dataframes/from_matteo/developability_pca_raw_habib/developability_query-4.csv")
library(reticulate)
np <- import("numpy")
setwd("~/developability_project/plotting_dataframes/from_matteo/developability_pca_raw_habib/")
query_4_matrix <- np$load("developability_query-4_metrics_seq.npy")
class(query_4_matrix)

#organising data query_4 (all antibodies)

query_4_data = query_4_data %>%
  select(rowid, starts_with("AbOrig")) %>% #discard sequence and unnecessary infro (for now)
  rename(originala_data_rowid = rowid) %>% #this is the original name of the antibody from the original data
  rowid_to_column() #add rowid to glue that matrix with the metadata as per Matteo's message " I didn't specify it so far but ofc the Abs in the corresponding npy & csv files are ordered in the same way" - 20221115 11:43 AM

query_4_data = query_4_data %>%
  rename(original_data_rowid = originala_data_rowid)


query_4_PCAs = as_tibble(query_4_matrix) %>%
  select(V1,V2) %>%
  rename(PCA1 = V1, PCA2 = V2) %>%
  rowid_to_column() 

query_4_plotting = query_4_PCAs %>%
  left_join(query_4_data) 


query_4_plotting = query_4_plotting %>%
  mutate(AbOrig_species = str_to_title(AbOrig_species))%>%
  rename(Species = AbOrig_species)%>%
  rename(Chain = AbOrig_chain) %>%
  mutate(Chain = str_to_title(Chain))

#this is panel C left
gscatter_sample_all = ggplot(query_4_plotting %>% sample_n(nrow(query_4_plotting)), # This strategy provides a way to shuffle the dataframe so no subset/species can hide the other (overlap over it). 
                             aes(x = PCA1, y = PCA2)) +
  geom_point(aes(fill = Chain, shape = Chain, colour = Chain), size = 0.5, alpha = 0.3) + 
  scale_colour_manual(values = c("#31a354", "#d9f0a3")) +
  scale_fill_manual(values = c("#31a354", "#d9f0a3")) +
  theme_bw() +
  xlab("PC1 (17.2%)") + ylab("PC2 (15.4%)")+
  scale_shape_manual(values=21:22) +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  theme(legend.position = "top",
    legend.justification = "left",
    legend.title = element_text(size=24),
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 16))


ggsave(gscatter_sample_all, file = "~/developability_project/plots/PCAs/proper_all_anative_antibodies_DP_embedings_scatter_chain.png", height = 8, width = 10)

#organising data query_2 (VH antibodies)

query_2_data = read_csv("~/developability_project/plotting_dataframes/from_matteo/developability_pca_raw_habib/developability_query-2.csv")
setwd("~/developability_project/plotting_dataframes/from_matteo/developability_pca_raw_habib/") #where the .npy frames are
query_2_matrix <- np$load("developability_query-2_metrics_seq.npy")

query_2_data = query_2_data %>%
  select(rowid, starts_with("AbOrig")) %>% #discard sequence and unnecessary info (for now)
  rename(original_data_rowid = rowid) %>% #this is the original name of the antibody from the original data
  rowid_to_column() 


query_2_PCs = as_tibble(query_2_matrix) %>%
  select(V1,V2) %>%
  rename(PC1 = V1, PC2 = V2) %>%
  rowid_to_column() 


query_2_plotting = query_2_PCs %>%
  left_join(query_2_data) %>%
  mutate(AbOrig_species = str_to_title(AbOrig_species)) %>%
  rename(Species = AbOrig_species) %>%
  rename(Chain = AbOrig_chain) %>%
  mutate(Chain = str_to_title(Chain)) 

#this is panel C top right

gscatter_all_q2 = ggplot(query_2_plotting  %>% sample_n(nrow(query_2_plotting)), aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Species, shape = Species, colour = Species), size = 0.5, alpha = 0.6) + 
  scale_colour_manual(values = c("#fdcdac", "#cbd5e8")) +
  scale_fill_manual(values = c("#fdcdac", "#cbd5e8")) +
  theme_bw() +
  xlab("PC1 (20.16%)") + ylab("PC2 (10.12%)")+
  scale_shape_manual(values=21:22) +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  theme(legend.position = "top",
    legend.justification = "left",
    legend.title = element_text(size=24),
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 20))

ggsave(gscatter_all_q2, file = "~/developability_project/plots/PCAs/proper_heavy_anative_antibodies_q3_DP_embedings_scatter_species.png", height = 8, width = 10)



# query-3 (VL antibodies)

query_3_data = read_csv("~/developability_project/plotting_dataframes/from_matteo/developability_pca_raw_habib/developability_query-3.csv")
setwd("~/developability_project/plotting_dataframes/from_matteo/developability_pca_raw_habib/")
query_3_matrix <- np$load("developability_query-3_metrics_seq.npy")

query_3_data = query_3_data %>%
  select(rowid, starts_with("AbOrig")) %>% #discard sequence and unnecessary infro (for now)
  rename(original_data_rowid = rowid) %>% #this is the original name of the antibody from the original data
  rowid_to_column() 


query_3_PCs = as_tibble(query_3_matrix) %>%
  select(V1,V2) %>%
  rename(PC1 = V1, PC2 = V2) %>%
  rowid_to_column() 


query_3_plotting = query_3_PCs %>%
  left_join(query_3_data) %>%
  mutate(AbOrig_species = str_to_title(AbOrig_species)) %>%
  rename(Species = AbOrig_species) %>%
  rename(Chain = AbOrig_chain) %>%
  mutate(Chain = str_to_title(Chain)) 

#this is panel C bottom right

gscatter_all_q3 = ggplot(query_3_plotting %>% sample_n(nrow(query_3_plotting)),
                         #sample_n(10),
                         aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Species, shape = Species, colour = Species), size = 0.5, alpha = 0.6) + 
  #geom_density_2d(aes(color = Chain)) +
  scale_colour_manual(values = c("#fdcdac", "#cbd5e8")) +
  scale_fill_manual(values = c("#fdcdac", "#cbd5e8")) +
  theme_bw() +
  xlab("PC1 (20.29%)") + ylab("PC2 (11.49%)")+
  scale_shape_manual(values=21:22) +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  theme(#legend.position=c(0.08,0.9),
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_text(size=24),
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 20))



ggsave(gscatter_all_q3, file = "~/developability_project/plots/PCAs/proper_light_anative_antibodies_q3_DP_embedings_scatter_species.png", height = 8, width = 10)
ggsave("~/developability_project/plots/PCAs/proper_light_anative_antibodies_q3_DP_embedings_scatter_species_legend.pdf", height = 8, width = 10)


#box/violin plots for Supp figure 4C-----

query_3_plotting_long = query_3_plotting %>% 
  select(PC1, PC2, Species,original_data_rowid,rowid) %>%
  pivot_longer(-c(Species,original_data_rowid,rowid), values_to = "PC", names_to = "componenet") %>%
  mutate(componenet = str_remove_all(componenet, "A"))

query_3_plotting_long_text = query_3_plotting_long %>%
  group_by(Species, componenet) %>%
  summarise(MD = median(PC))


gbox_violin_q3 = ggplot(query_3_plotting_long, aes(x = Species, y = PC)) +
  facet_grid(cols = vars(componenet)) +
  geom_violin(aes(fill = Species), alpha = 0.5) +
  geom_boxplot(aes(fill = Species), outlier.shape = "", width = 0.35) +
  scale_fill_manual(values = c("#fdcdac", "#cbd5e8")) +
  theme_bw() +
  geom_text(data = query_3_plotting_long_text, aes(x = Species, y = MD+0.65, label = round(MD, digits = 1)), size = 3.75) +
  xlab("") + ylab("PC value") +
  theme(strip.background = element_rect(fill = "white"), 
        strip.text = element_text(size = 16), 
        axis.title = element_text(size = 16),
        axis.text.x = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        axis.ticks.x = element_blank())

ggsave(gbox_violin_q3, file = "~/developability_project/plots/PCAs/light_anative_antibodies_DP_embedings_boxviolin_species.png", height = 6, width = 4)
ggsave(gbox_violin_q3, file = "~/developability_project/plots/PCAs/light_anative_antibodies_DP_embedings_boxviolin_species.pdf", height = 6, width = 4)



query_2_plotting_long = query_2_plotting %>% 
  select(PC1, PC2, Species,original_data_rowid,rowid) %>%
  pivot_longer(-c(Species,original_data_rowid,rowid), values_to = "PC", names_to = "componenet") %>%
  mutate(componenet = str_remove_all(componenet, "A"))


query_2_plotting_long_text = query_2_plotting_long %>%
  group_by(Species, componenet) %>%
  summarise(MD = median(PC))

gbox_violin_q2 = ggplot(query_2_plotting_long, aes(x = Species, y = PC)) +
  facet_grid(cols = vars(componenet)) +
  geom_violin(aes(fill = Species), alpha = 0.5) +
  geom_boxplot(aes(fill = Species), outlier.shape = "", width = 0.35) +
  scale_fill_manual(values = c("#fdcdac", "#cbd5e8")) +
  theme_bw() +
  geom_text(data = query_2_plotting_long_text, aes(x = Species, y = MD+0.65, label = round(MD, digits = 1)), size = 3.75) +
  xlab("") + ylab("PC value") +
  theme(strip.background = element_rect(fill = "white"), 
        strip.text = element_text(size = 16), 
        axis.title = element_text(size = 16),
        axis.text.x = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        axis.ticks.x = element_blank())

ggsave(gbox_violin_q2, file = "~/developability_project/plots/PCAs/heavy_anative_antibodies_DP_embedings_boxviolin_species.png", height = 6, width = 4)
ggsave(gbox_violin_q2, file = "~/developability_project/plots/PCAs/heavy_anative_antibodies_DP_embedings_boxviolin_species.pdf", height = 6, width = 4)



query_4_plotting_long = query_4_plotting %>% 
  select(PCA1, PCA2, Chain,original_data_rowid,rowid) %>%
  pivot_longer(-c(Chain,original_data_rowid,rowid), values_to = "PC", names_to = "componenet") %>%
  mutate(componenet = str_remove_all(componenet, "A"))

query_4_plotting_long_text = query_4_plotting_long %>%
  group_by(Chain, componenet) %>%
  summarise(MD = median(PC))

gbox_violin_all = ggplot(query_4_plotting_long, aes(x = Chain, y = PC)) +
  facet_grid(cols = vars(componenet)) +
  geom_violin(aes(fill = Chain), alpha = 0.5) +
  geom_boxplot(aes(fill = Chain), outlier.shape = "", width = 0.35) +
  scale_fill_manual(values = c("#31a354", "#d9f0a3")) +
  theme_bw() +
  geom_text(data = query_4_plotting_long_text, aes(x = Chain, y = MD+0.65, label = round(MD, digits = 1)), size = 3.75) +
  xlab("") + ylab("PC value") +
  theme(strip.background = element_rect(fill = "white"), 
        strip.text = element_text(size = 16), 
        axis.title = element_text(size = 16),
        axis.text.x = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        axis.ticks.x = element_blank()) 

ggsave(gbox_violin_all, file = "~/developability_project/plots/PCAs/all_anative_antibodies_DP_embedings_boxviolin_chain.png", height = 6, width = 4)
ggsave(gbox_violin_all, file = "~/developability_project/plots/PCAs/all_anative_antibodies_DP_embedings_boxviolin_chain.pdf", height = 6, width = 4)


