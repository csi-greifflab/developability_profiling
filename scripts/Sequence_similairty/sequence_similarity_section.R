#Author: Habib Bashour

library(seriation)
library(janitor)
library(reticulate)
library(reshape2)

read.csv("AbChain_native_developability_final.csv") #from data/native directory 
read.csv("AbStruc_native_developability_final.csv") #from data/native directory 



#try to take a random sample from all the antibodies in the isotype/species instead of just the ones with complete immunogenecity reads

#first create new directory structure------ 

setwd("~/developability_project/environments_DataObjects")
dir.create(file.path("dpc_matrices_same_vgene"), recursive = TRUE) #it works
setwd("~/developability_project/environments_DataObjects/dpc_matrices_same_vgene/")
dir.create(file.path("sequence"), recursive = TRUE) #it works
dir.create(file.path("structure"), recursive = TRUE) #it works

setwd("~/developability_project/environments_DataObjects/dpc_matrices_same_vgene/sequence/")
batch100_dirs = c(as.character(seq(1,100,1)))
for(d in 1:length(batch100_dirs)){
  dir.create(file.path(batch100_dirs[d]), recursive = TRUE)
}

setwd("~/developability_project/environments_DataObjects/dpc_matrices_same_vgene/structure/")
batch100_dirs = c(as.character(seq(1,100,1)))
for(d in 1:length(batch100_dirs)){
  dir.create(file.path(batch100_dirs[d]), recursive = TRUE)
}



setwd("~/developability_project/environments_DataObjects/dpc_matrices_same_vgene/")
dir.create(file.path("batch_rowids"), recursive = TRUE) #it works
setwd("~/developability_project/environments_DataObjects/dpc_matrices_same_vgene/batch_rowids/")
for(d in 1:length(batch100_dirs)){
  dir.create(file.path(batch100_dirs[d]), recursive = TRUE)
}




#scale the three main dataframes------ 
scale_manual_fun_na = function(x){
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

AbChain_native_developability_scaled %>%sample_n(100) %>% View
AbChain_native_developability_scaled = AbChain_native_developability_final %>% #200K
  select(-c( database,starts_with("aaSeq") ,v_gene ,j_gene)) %>%
  select(-c( AbChain_immunopeptide_end_region ,AbChain_immunopeptide_start_region)) %>% #remove categorical variables
  group_by(chain,species) %>% #customise the scaling/normalisation by chain
  mutate_at(scale_manual_fun_na, .vars = vars(-c(chain,species, isotype, rowid)))

AbStruc_native_developability_scaled = AbStruc_native_developability_final %>%
  select(-c(database,starts_with("aaSeq") ,v_gene ,j_gene)) %>%
  select(-c(AbStruc_weak_hbonds, AbStruc_pi_helices)) %>%
  select(where(~n_distinct(.) > 1)) %>% #drop columns (parameters) with same value across all sequences
  group_by(chain,species) %>% #customise the scaling/normalisation by chain 
  mutate_at(scale_manual_fun_na, .vars = vars(-c(chain,species, isotype, rowid)))



save(AbChain_native_developability_scaled, file = "~/developability_project/environments_DataObjects/AbChain_native_developability_scaled.RData")
save(AbStruc_native_developability_scaled, file = "~/developability_project/environments_DataObjects/AbStruc_native_developability_scaled.RData")

load("~/developability_project/environments_DataObjects/AbChain_native_developability_scaled.RData")
load("~/developability_project/environments_DataObjects/AbStruc_native_developability_scaled.RData")

#load scaled dataframes and vgene representation info------ 

retrieve_vgene = AbChain_native_developability_final %>%
  select(rowid,v_gene)



#Run a loop that can calculate dpc for thr three levels and sequence similarity for the same batch of sequences (100 abs X 100 times)------ 

SpecieS = c("human", "human" ,"human", "human", "human", "human" ,"human" ,"mouse", "mouse" ,"mouse" ,"mouse")
isotypes = c("IgG", "IgA" ,"IgM" ,"IgD", "IgE", "IgK", "IgL" ,"IgG", "IgM", "IgK" ,"IgL")

#run the following line before the loop to empty the vector
used_rowids_vgenes = c()

for (h in 1:length(isotypes)){
  for (s in 1:length(batch100_dirs)){
    taget_vgene_family = v_gene_family_distibution_native %>%
      drop_na() %>%
      filter(!str_detect(v_gene_family, "Musspr")) %>% 
      filter(n>5000) %>%
      #view() %>%
      filter(species == SpecieS[h],
             isotype == isotypes[h]) %>% 
      pull(v_gene_family) %>%
      sample(1)
    target_data_seq = AbChain_native_developability_scaled %>%
      filter(isotype == isotypes[h]) %>%
      filter(species == SpecieS[h]) %>%
      mutate(v_gene_family = str_remove_all(v_gene, "\\-.*")) %>%
      ungroup() %>%
      filter(!rowid %in% used_rowids_vgenes) %>%
      filter(v_gene_family == taget_vgene_family) %>%
      select(-c(isotype,species,chain,v_gene,v_gene_family)) %>%
      sample_n(100) %>% 
      arrange(rowid) %>% #it is of extreme importance to arrange so I avoid the lose of data when binding melted matrices (lost half a day here).
      t() %>% #flip the rows to columns
      as_tibble() %>% 
      row_to_names(row_number = 1) %>%
      mutate_all(~ as.numeric(.))
    used_in_target_data = as.double(colnames(target_data_seq)) 
    used_rowids_vgenes = c(used_in_target_data, used_rowids_vgenes)
    AbChain_100_seq_seq_pearson_cor_matrix_batch =  round(cor(target_data_seq, method = "pearson", use = "pairwise.complete.obs"),2) 
    save(AbChain_100_seq_seq_pearson_cor_matrix_batch, file = paste0("~/developability_project/environments_DataObjects/dpc_matrices_same_vgene/sequence/", batch100_dirs[s],"/AbChain_",isotypes[h],"_",SpecieS[h],"_100_seq_seq_pearson_dpc_matrix_batch_",batch100_dirs[s],".RData"))
    saveRDS(used_in_target_data, file = paste0("~/developability_project/environments_DataObjects/dpc_matrices_same_vgene/batch_rowids/",batch100_dirs[s],"/batch_",batch100_dirs[s],"_",isotypes[h],"_",SpecieS[h],"_rowids.rds"))
    
    target_data_struc = AbStruc_native_developability_scaled %>%
      filter(rowid %in% used_in_target_data) %>%
      arrange(rowid) %>%
      select(where(~!any(is.na(.)))) %>% #drop columns with at least one NaN
      ungroup() %>%
      select(-c(isotype,species,chain)) %>% 
      t() %>% #flip the rows to columns
      as_tibble() %>% 
      row_to_names(row_number = 1) 
    
    AbStruc_100_seq_seq_pearson_cor_matrix_batch =  round(cor(target_data_struc, method = "pearson", use = "pairwise.complete.obs"),2) 
    save(AbStruc_100_seq_seq_pearson_cor_matrix_batch, file = paste0("~/developability_project/environments_DataObjects/dpc_matrices_same_vgene/structure/", batch100_dirs[s],"/AbStruc_",isotypes[h],"_",SpecieS[h],"_100_seq_seq_pearson_dpc_matrix_batch_",batch100_dirs[s],".RData"))
    print(paste0("done all DPC batch ",as.character(batch100_dirs[s]), " ",isotypes[h], " ", SpecieS[h]))
    
  }
}


#save all sampled rowids as a group 
saveRDS(used_rowids_vgenes, file = paste0("~/developability_project/environments_DataObjects/dpc_matrices_same_vgene/batch_rowids/all_sampled_rowids_vgenes.rds"))



#now sequence similarity matrices------


setwd("~/developability_project/environments_DataObjects/seq_seq_similarity_matrices/same_vgene_familiy_100batches/")
for(d in 1:length(batch100_dirs)){
  dir.create(file.path(batch100_dirs[d]), recursive = TRUE)
}


stringsimilarityFast_habib = function(word1, word2){
  1 - (stringdist(word1, word2, method = "lv")/max(nchar(word1), nchar(word2)))
}

StringSimilarityMat_habib <- function(seq_rowid_data){
  string_vector = seq_rowid_data$aaSeqAbChain
  names_vector = seq_rowid_data$rowid
  m = diag(0, length(string_vector))
  sapply(1:(length(string_vector)-1), function(i){ #apply this function i number where i = dim -1 (exclude the diagonal)
    m[,i] <<- c(rep(0,i),lapply(string_vector[(i+1):length(string_vector)],stringsimilarityFast_habib,string_vector[i])) %>% unlist() 
  }) #populate only columns
  m = m+t(m) #produce a symetric matrix
  dimnames(m) = list(names_vector,names_vector)
  diag(m) <- 1 # set the diagonal to 1 (similarity of a sequence with itself)
  return(m)
}

for(B in 1:length(batch100_dirs)){
  setwd(paste0("~/developability_project/environments_DataObjects/dpc_matrices_same_vgene/batch_rowids/",batch100_dirs[B],"/")) #we can take the rowids from any correlation matrix
  list_rds = dir(pattern = ".rds") #we only need one of them to retrieve the rowids
  for (M in 1:length(list_rds)){
    batch_rowids = readRDS(paste0("~/developability_project/environments_DataObjects/dpc_matrices_same_vgene/batch_rowids/",batch100_dirs[B],"/",list_rds[M]))
    category = str_remove(list_rds[M], paste0("batch_", "[:digit:]+", "_"))
    category = str_remove(category, "_rowids.rds")
    
    file_name = str_replace(list_rds[M], "rowids.rds", "seq_sim_matrix")
    
    sample_seqs_rowids = AbChain_native_developability_final %>%
      filter(rowid %in% batch_rowids) %>% 
      arrange(rowid) %>%
      select(rowid, aaSeqAbChain)
    #produce the similarity matrix using the custom function StringSimilarityMat_habib
    similarity_matrix = StringSimilarityMat_habib(sample_seqs_rowids)
    saveRDS(similarity_matrix, file = paste0("~/developability_project/environments_DataObjects/seq_seq_similarity_matrices/same_vgene_familiy_100batches/",batch100_dirs[B],"/",file_name,".rds"))
    print(paste0(file_name, " done!"))
  }
}




#produce heatmaps in Figure 5A and supp figures 6,7,8-----

directories = c("sequence", "structure")

random_batches =  round(runif(10, 1, 100), digits=0)
for (RB in length(random_batches)){
  for(D in 1:length(directories)){
    setwd(paste0("~/developability_project/environments_DataObjects/dpc_matrices_same_vgene/", directories[D], "/",random_batches[RB],"/"))
    list_RData_cor = dir(pattern = ".RData")
    setwd(paste0("~/developability_project/environments_DataObjects/seq_seq_similarity_matrices/same_vgene_familiy_100batches/",random_batches[RB],"/"))
    list_RData_sim = dir(pattern = "rds")
    
    #check if the names are ordered similarly --> check is done --> similar order
    for (M in 1:length(list_RData_cor)){
      #load the matrix
      cor_matrix = loadRData(paste0("~/developability_project/environments_DataObjects/dpc_matrices_same_vgene/", directories[D], "/",random_batches[RB],"/",list_RData_cor[M]))
      sim_matrix = readRDS(paste0("~/developability_project/environments_DataObjects/seq_seq_similarity_matrices/same_vgene_familiy_100batches/",random_batches[RB],"/",list_RData_sim[M]))
      plot_name = str_replace_all(list_RData_sim[M], "matrix.rds", "vs_DPC")
      #order the seq_seq correlation matrix using hirachial clustering; complete link
      cor_matrix_oder = seriate(dist(cor_matrix), method = "HC_complete")
      #vgene annotation is not needed anymore because we have all the antibodies on the same batch from same gene family
      
      #plotting------
      pdf(paste0("~/developability_project/plots/seq_seq_100batch_correlation_Vs_similarity/",directories[D],"/",plot_name, ".pdf"), height = 8,width = 14)
      diag(cor_matrix) <- NA
      ht1 = Heatmap(cor_matrix, name = "DPC (Pearson)",
                    col = circlize::colorRamp2(c(-1,0,1), c("#b2182b", "#ffffff", "#4d4d4d")),
                    column_order = get_order(cor_matrix_oder),
                    row_order = get_order(cor_matrix_oder),
                    show_row_dend = F, show_column_dend = F,
                    show_heatmap_legend = T,
                    #use_raster = TRUE,
                    #raster_resize = T, raster_device = "png",
                    #column_title = "developability profile correlation",
                    #right_annotation = c(vgene_annotation_right),
                    na_col = "white",
                    row_names_gp = gpar(fontsize = 0.1),column_names_gp = gpar(fontsize = 0.1))
      
      diag(sim_matrix) <- NA
      
      ht2 = Heatmap(sim_matrix, name = "sequence similarity score",
                    #col = circlize::colorRamp2(c(min(sim_matrix, na.rm = T),max(sim_matrix, na.rm = T)), c("white", "#238b45")),
                    col = circlize::colorRamp2(c(quantile(sim_matrix, prob =.05, na.rm = T),max(sim_matrix, na.rm = T)), c("white", "#238b45")), #only batch 57 was done like this
                    
                    column_order = get_order(cor_matrix_oder),
                    row_order = get_order(cor_matrix_oder),
                    show_row_dend = F, show_column_dend = F,
                    show_heatmap_legend = T,
                    na_col = "white",
                    #raster_resize = T, raster_device = "png",
                    #column_title = "sequence similarity",
                    #right_annotation = c(vgene_annotation_right),
                    row_names_gp = gpar(fontsize = 0),column_names_gp = gpar(fontsize = 0))
      
      ht_list = ht1 + ht2
      draw(ht_list)
      dev.off()
      print(paste0("done batch ", plot_name,  " !"))
    }
  }
  
}




#produce DPC and sequence similarity bar plot in Figure 5B for 100 batches -------


batches100_cor_sim_correlation = tibble(category= character(), kind = character(), average_seq_similaritiy = as.double(),
                                        batch = as.character(), pearson = numeric(), abs_pearson = numeric(),
                                        spearman = numeric(), abs_spearman = numeric())




batch100_dirs = c(as.character(seq(1,100,1)))

View(batches100_cor_sim_correlation)
for(D in 1:length(directories)){
  for(B in 1:length(batch100_dirs)){
    setwd(paste0("~/developability_project/environments_DataObjects/dpc_matrices_same_vgene/", directories[D], "/",batch100_dirs[B],"/"))
    list_RData_cor = dir(pattern = "dpc_matrix_batch_.+.RData") #
    setwd(paste0("~/developability_project/environments_DataObjects/seq_seq_similarity_matrices/same_vgene_familiy_100batches/",batch100_dirs[B],"/"))
    list_RData_sim = dir(pattern = "sim_matrix.rds")
    categories = str_remove_all(list_RData_cor, "AbChain_|AbStruc_|AbChStruc_") %>%
      str_remove_all("_100_seq_seq_.+.RData")
    
    for(M in 1:length(list_RData_cor)){
      cor_matrix = loadRData(paste0("~/developability_project/environments_DataObjects/dpc_matrices_same_vgene/",directories[D],"/",batch100_dirs[B],"/",list_RData_cor[M]))
      sim_matrix = readRDS(paste0("~/developability_project/environments_DataObjects/seq_seq_similarity_matrices/same_vgene_familiy_100batches/",batch100_dirs[B],"/",list_RData_sim[M]))
      
      cor_matrix[lower.tri(cor_matrix, diag = TRUE)] <- NA
      sim_matrix[lower.tri(sim_matrix, diag = TRUE)] <- NA
      mean_sequence_similarity_score = mean(sim_matrix, na.rm = T)
      cor_df = melt(cor_matrix) %>%
        as_tibble() %>%
        drop_na() %>%
        mutate(category = categories[M]) %>%
        mutate(kind = directories[D]) %>%
        mutate(batch = batch100_dirs[B]) %>%
        rename(pearson = value)
      sim_df = melt(sim_matrix) %>%
        as_tibble() %>%
        drop_na() %>%
        mutate(category = categories[M]) %>%
        mutate(kind = directories[D]) %>%
        mutate(batch = batch100_dirs[B]) %>%
        rename(sim_score = value)
      
      cor_sim_df = left_join(cor_df, sim_df) %>%
        select(-c(Var1,Var2))
      cor_sim_coR_pearson = cor(cor_sim_df$pearson, cor_sim_df$sim_score)
      abs_cor_sim_coR_pearson= cor(abs(cor_sim_df$pearson), abs(cor_sim_df$sim_score))
      cor_sim_coR_spearman = cor(cor_sim_df$pearson, cor_sim_df$sim_score, method = "spearman")
      abs_cor_sim_coR_spearman = cor(abs(cor_sim_df$pearson), abs(cor_sim_df$sim_score), method = "spearman")
      
      
      batches100_cor_sim_correlation = add_row(batches100_cor_sim_correlation, category = categories[M], kind = directories[D],
                                               batch = batch100_dirs[B], average_seq_similaritiy = mean_sequence_similarity_score,
                                               pearson = cor_sim_coR_pearson, abs_pearson = abs_cor_sim_coR_pearson,
                                               spearman = cor_sim_coR_spearman, abs_spearman = abs_cor_sim_coR_spearman)
      
      #cor_sim_df_final = rbind(cor_sim_df_final, cor_sim_df)
    }
  }
  
}

pearson_sum_100batch = batches100_cor_sim_correlation %>%
  filter(measure %in% c("average_seq_similaritiy", "pearson")) %>%
  group_by(isotype, species, kind, measure) %>%
  summarise(mean_value = mean(value), SD_value = sd(value))


isotype_plotting_order = c("IgD","IgM","IgG","IgA","IgE","IgK","IgL")

pearson_sum_100batch$kind = factor(pearson_sum_100batch$kind, levels = c("sequence", "structure", "sequence_and_structure"))
pearson_sum_100batch$isotype = factor(pearson_sum_100batch$isotype, levels = c("IgD","IgM","IgG","IgA","IgE","IgK","IgL"))
pearson_sum_100batch$measure = factor(pearson_sum_100batch$measure, levels = c("pearson", "average_seq_similaritiy"))
pearson_sum_100batch = pearson_sum_100batch %>%
  mutate(species = str_to_title(species),
         kind = str_to_title(kind))

ggplot(pearson_sum_100batch %>% filter(kind!= "Sequence_and_structure"), aes(x = isotype, y = mean_value, fill = measure)) +
  facet_grid(cols = vars(species), rows = vars(kind), scales = "free_x", space='free_x')+
  geom_bar(position = position_dodge(), stat = "identity", alpha = 0.75,width=0.8) +
  geom_errorbar(aes(ymin = mean_value - SD_value, ymax = mean_value + SD_value), position=position_dodge(width=0.9), width = 0.3, colour="black")+
  theme_bw() +
  geom_text(aes(label = round(mean_value, digits = 1), y = - 0.1, group = measure, colour = measure), position = position_dodge(width = 0.9), size = 3.3, show.legend = F) +
  #geom_text(aes(label = round(mean_value, digits = 1), y = -0.1, group = measure, colour = measure), position = position_dodge(width = 0.9), size = 3.3, show.legend = F) +
  geom_point(data =  batches100_cor_sim_correlation %>% filter(measure %in% c("average_seq_similaritiy", "pearson")) %>% filter(kind!= "Sequence_and_structure"),
             aes(fill=measure ,x = isotype, y = value),
             position = position_jitterdodge(jitter.height = 0.07, dodge.width = 1),
             pch = 21, size = 0.8, colour = "black", stroke = 0.3, alpha = 0.95,
             show.legend = FALSE) +
  xlab("")+ ylab("Pearson correlation and average sequence similarity") +
  theme(strip.background = element_rect(fill = "white"),
        legend.position="top", legend.justification='left',
        strip.text = element_text(size = 14), plot.title =  element_text(size = 12),
        plot.caption = element_text(face = "italic", hjust = 0),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  #ggtitle("seq_seq developability landscape sequence similarity correlation Vs sequence similarity", subtitle = "10 batches; 1000 seqs each") +
  scale_fill_manual(name = "Metric", labels = c("Pearson correlation", "Sequence similarity"),
                    values = c("#d8b365","#4e9d60")) +
  scale_colour_manual(values=c("#d8b365","#4e9d60"))


ggsave2("~/developability_project/plots/seq_seq_100batch_correlation_Vs_similarity/cor_sim_pcorrelation_100batches_ready_panel.png", width = 10, height = 6)
ggsave2("~/developability_project/plots/seq_seq_100batch_correlation_Vs_similarity/cor_sim_pcorrelation_100batches_ready_panel.pdf", width = 10, height = 6)





#for Supp Figure 9 panel A------
batch =  round(runif(1, 1, 100), digits=0) #choose one number between 1-100

example_plotting_empty= left_join(get_melted_matrix(batch_seqsim_matrix, "seq_sim"), get_melted_matrix(batch_sequence_matrix, "seq_dpc")) %>%
  left_join(get_melted_matrix(batch_structure_matrix, "struc_dpc")) %>%
  pivot_longer(-c(rowid1 , rowid2 ,seq_sim), values_to = "DPC", names_to = "level") %>%
  mutate(level = str_replace_all(level, "seq_dpc", "sequence")) %>%
  mutate(level = str_replace_all(level, "struc_dpc", "structure"))%>%
  mutate(level = str_replace_all(level, "seqstructure", "sequence and structure")) %>%
  mutate(category == "") %>%
  filter(category == "habib")

#useful loop to plot seqsim vs dpc values (the seqsim vs dpc loop)-------
#batch = 83
for (i in 1:11){
  setwd(paste0("~/developability_project/environments_DataObjects/dpc_matrices_same_vgene/sequence/",batch,"/"))
  list_RData_seuqence = dir(pattern = ".RData")
  setwd(paste0("~/developability_project/environments_DataObjects/dpc_matrices_same_vgene/structure/",batch,"/"))
  list_RData_structure = dir(pattern = ".RData")
  setwd(paste0("~/developability_project/environments_DataObjects/seq_seq_similarity_matrices/same_vgene_familiy_100batches/",batch,"/"))
  list_RData_seqsim = dir(pattern = ".rds")
  
  batch_sequence_matrix = loadRData(paste0("~/developability_project/environments_DataObjects/dpc_matrices_same_vgene/sequence/",batch,"/",list_RData_seuqence[i]))
  batch_structure_matrix = loadRData(paste0("~/developability_project/environments_DataObjects/dpc_matrices_same_vgene/structure/",batch,"/",list_RData_structure[i]))
  batch_seqsim_matrix = readRDS(paste0("~/developability_project/environments_DataObjects/seq_seq_similarity_matrices/same_vgene_familiy_100batches/",batch,"/",list_RData_seqsim[i]))
  
  category = str_remove_all(list_RData_seuqence[i],"AbChain_")
  category = str_remove_all(category, "_100_seq_seq_.*")
  
  file_plotting = left_join(get_melted_matrix(batch_seqsim_matrix, "seq_sim"), get_melted_matrix(batch_sequence_matrix, "seq_dpc")) %>%
    left_join(get_melted_matrix(batch_structure_matrix, "struc_dpc")) %>%
    left_join(get_melted_matrix(batch_both_matrix, "seqstruc_dpc")) %>%
    pivot_longer(-c(rowid1 , rowid2 ,seq_sim), values_to = "DPC", names_to = "level") %>%
    mutate(level = str_replace_all(level, "seq_dpc", "sequence")) %>%
    mutate(level = str_replace_all(level, "struc_dpc", "structure"))%>%
    mutate(category = category)
  
  example_plotting_empty = rbind(file_plotting, example_plotting_empty)
  
}


IgM_83_seqsim_vs_dpc = example_plotting_empty %>%
  filter(category == "IgM_human") %>%
  filter(level != "sequence and structure") %>%
  mutate(level = str_to_title(level))



scatter_and_density(IgM_83_seqsim_vs_dpc %>% filter(level == "Sequence"), x = "seq_sim", y = "DPC", title_input = "Sequence")
ggsave("~/developability_project/plots/seq_seq_100batch_correlation_Vs_similarity/IgM_83_seqsim_vs_seqdpc.png", height = 6, width = 6)
ggsave("~/developability_project/plots/seq_seq_100batch_correlation_Vs_similarity/IgM_83_seqsim_vs_seqdpc.pdf", height = 6, width = 6)

scatter_and_density(IgM_83_seqsim_vs_dpc %>% filter(level == "Structure"), x = "seq_sim", y = "DPC", title_input = "Structure")
ggsave("~/developability_project/plots/seq_seq_100batch_correlation_Vs_similarity/IgM_83_seqsim_vs_strucdpc.png", height = 6, width = 6)
ggsave("~/developability_project/plots/seq_seq_100batch_correlation_Vs_similarity/IgM_83_seqsim_vs_strucdpc.pdf", height = 6, width = 6)



#produce PCA plot for sequence similarity clusters in Figure 5C---------


query_1_data = read_csv("~/developability_project/plotting_dataframes/from_matteo/developability_query-1/developability_query-1.csv")
np <- import("numpy")
setwd("~/developability_project/plotting_dataframes/from_matteo/developability_query-1/")

query_1_matrix_full_DPL <- np$load("developability_query-1_metrics_seq.npy")


query_1_full_DPL_PCAs = as_tibble(query_1_matrix_full_DPL) %>%
  select(V1,V2) %>%
  rename(PC1 = V1, PC2 = V2) %>%
  rowid_to_column() 

query_1_full_DPL_cluster_plotting = query_1_full_DPL_PCAs %>%
  left_join(query_1_data) 

query_1_full_DPL_cluster_plotting = query_1_full_DPL_cluster_plotting %>% 
  mutate(AbChain_big_clusters_75_names = case_when(AbChain_big_clusters_75 == 1182371 ~ "1",
                                                   AbChain_big_clusters_75 == 1183318 ~ "2",
                                                   AbChain_big_clusters_75 == 1182400 ~ "3",
                                                   AbChain_big_clusters_75 == 1182384 ~ "4",
                                                   AbChain_big_clusters_75 == 1182672 ~ "5",
                                                   AbChain_big_clusters_75 == 1183315 ~ "6",
                                                   AbChain_big_clusters_75 == 1267619 ~ "7",
                                                   TRUE ~ "non-clustered")) %>%
  rename(`Cluster ID` = AbChain_big_clusters_75_names)


x = query_1_full_DPL_cluster_plotting %>% filter(`Cluster ID` == "non-clustered")
y =   query_1_full_DPL_cluster_plotting %>% filter(`Cluster ID` != "non-clustered")

clusters_full_DPL_plot = ggplot(rbind(x,y), 
                                #clusters_full_DPL_plot = ggplot(query_1_full_DPL_cluster_plotting %>% sample_n(nrow(query_1_full_DPL_cluster_plotting)), 
                                #%>% sample_n(10), only for the legend plotting
                                aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = `Cluster ID`, colour = `Cluster ID`), size = 0.3, alpha = 0.3) + 
  #geom_density_2d(aes(color = Chain)) +
  #scale_color_brewer(palette = "Accent")+
  #scale_fill_brewer(palette = "Accent")+
  scale_colour_manual(values = c("#7fc97f", "#beaed4", "#fdc086", "#ffff99","#386cb0","#f0027f","#bf5b17", "#666666")) +
  scale_fill_manual(values = c("#7fc97f", "#beaed4", "#fdc086", "#ffff99","#386cb0","#f0027f","#bf5b17", "snow")) +
  theme_bw() +
  xlab("PC1 (13.8%)") + ylab("PC2 (10.3%)")+ #provided by Matteo via Slack
  #scale_shape_manual(values=21:22) +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  theme(#legend.position=c(0.08,0.9),
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_text(size=24),
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 16))  




ggsave(clusters_full_DPL_plot, file = "~/developability_project/plots/PCAs/sequence_sim_clusters_arranged_full_DPL_embedding.png", height = 6, width = 10)  



#plot antibody count per sequence similarity cluster (Supp Figure 9B) -------

query_1_full_DPL_cluster_plotting %>%
  group_by(`Cluster ID`) %>%
  summarise(n = n()) %>%
  filter( `Cluster ID` != "non-clustered") %>%
  ggplot(aes(x = `Cluster ID`, y = n)) +
  geom_col(aes(fill = `Cluster ID`, colour = `Cluster ID`), alpha = 0.3) +
  scale_colour_manual(values = c("#7fc97f", "#beaed4", "#fdc086", "#ffff99","#386cb0","#f0027f","#bf5b17")) +
  scale_fill_manual(values = c("#7fc97f", "#beaed4", "#fdc086", "#ffff99","#386cb0","#f0027f","#bf5b17")) +
  ylab("Antibody count") + xlab("Group ID") + 
  theme_bw()+
  guides(color = guide_legend(override.aes = list(size = 10))) +
  theme(#legend.position=c(0.08,0.9),
    text=element_text(family="Helvetica"),
    legend.position = "none",
    #legend.justification = "left",
    legend.title = element_text(size=24),
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 16)) +
  geom_text(aes(x = `Cluster ID`, y = n+1000, label = n))


ggsave("~/developability_project/plots/PCAs/sequence_sim_clusters_count.png", height = 7, width = 8)  














# hexplot for Supp Figure 9C (Ed/LD correlation analysis)------


#loading numpy module
np <- import("numpy")
setwd("~/developability_project/plotting_dataframes/from_matteo/LD_ED_analysis/uncropped_slack_april_2023/complete-mwds_igm_iso-vgene-cluster-5k/")
ed_ditance_matrix_IgM <- np$load("sample_embs_dists.npy")
class(ed_ditance_matrix_IgM)
dim(ed_ditance_matrix_IgM)
ed_ditance_matrix_IgM[1:10,1:10]
#convert the symmetrical matrix to a df after removing the diagonal and the bottom triangle
ed_ditance_matrix_IgM[lower.tri(ed_ditance_matrix_IgM, diag = TRUE)] <- NA

ed_ditance_df_IgM = ed_ditance_matrix_IgM %>%
  melt() %>%
  drop_na() %>%
  dplyr::rename(ed = value)
rm(ed_ditance_matrix_IgM)  

ld_ditance_matrix_IgM <- np$load("sample_seqs_dists.npy")
ld_ditance_matrix_IgM[lower.tri(ld_ditance_matrix_IgM, diag = TRUE)] <- NA

ld_ditance_df_IgM = ld_ditance_matrix_IgM %>%
  melt() %>%
  drop_na() %>%
  dplyr::rename(ld = value)
rm(ld_ditance_matrix_IgM)  

ld_ditance_df_IgM = ld_ditance_df_IgM %>%
  as_tibble()

joined_df_IgM = ed_ditance_df_IgM %>%
  as_tibble() %>%
  right_join(ld_ditance_df_IgM)



joined_df_IgM %>%
  dplyr::rename(normalised_sequence_difference = ld) %>%
  mutate(sequence_similarity = 1- normalised_sequence_difference) %>% 
  filter(sequence_similarity != 1) %>%
  filter(euclidean_category == "close_proximity") %>%
  ggplot(aes(x = normalised_sequence_difference, y = ed)) +
  geom_hex(bins = 30, colour = "white", aes(fill = after_stat(count)))+
  stat_cor(aes(label = after_stat(r.label))) +
  theme_bw() +
  xlab("Sequence similarity") + 
  ylab("Euclidean distance") +
  theme(strip.background = element_rect(fill = "white"),
        text=element_text(family="Helvetica"),
        legend.text = element_text(size = 18),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        legend.title = element_text(size=22)) +
  scale_fill_viridis_c(name="Measurements",
                       breaks = c(50000, 100000, 150000), labels = c("50K", "100K", "150K")) 

ggsave("~/developability_project/plots/PCAs/ed_ld_close_proximity_hexed_uncropped.png",height = 7, width = 10) 
ggsave("~/developability_project/plots/PCAs/ed_ld_close_proximity_hexed_uncropped.pdf",height = 7, width = 10) 





