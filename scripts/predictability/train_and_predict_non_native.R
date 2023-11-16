library(tidyverse)
library(ggpubr)
library(data.table)
library(glue)
library(reticulate)

scale_manual_fun <- function(x) {
  (x - mean(x)) / sd(x)
}

scale_data = function(df) {
  df %>%
    mutate_at(scale_manual_fun, .vars = vars(starts_with(c("Ab"))))
}

Sys.setenv(RETICULATE_PYTHON = "python_env/bin/python") #Set python binary path as environment variable
np <- import("numpy")

# Import MWDS sets
AbStruc_mwds = readRDS('data/AbStruc_mwds.rds')
AbChain_mwds = readRDS('data/AbChain_mwds.rds')
AbChStruc_mwds = readRDS('data/AbChStruc_mwds.rds')

# Import DP data
AbChStruc_all = np$load("data/human-heavy_metrics_normed_updated.npy") %>%
  as_tibble()


# Import column names
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
  'AbStruc_pi_helices',
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
  'AbStruc_weak_hbonds',
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

# Add column names and index column
colnames(AbChStruc_all) = params
AbChStruc_all = AbChStruc_all %>%
  mutate(rowid = 1:nrow(.)) %>% 
  select(-c(AbStruc_pi_helices, AbStruc_steric_clashes, AbStruc_weak_hbonds))
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
  #'AbStruc_pi_helices',
  'AbStruc_total_interactions',
  'AbStruc_metsulphur_pi_interactions',
  'AbStruc_beta_bends',
  'AbStruc_beta_turns',
  'AbStruc_aromatic_interactions',
  'AbStruc_atom_atom_sum',
  'AbStruc_atom_plane_sum',
  'AbStruc_carbonyl_interactions',
  #'AbStruc_steric_clashes',
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
  #'AbStruc_weak_hbonds',
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
# Import test sets 

test_idxs = np$load("data/developability-test_idxs.npy")
native_data = AbChStruc_all %>%
  filter(rowid %in% test_idxs) %>% 
  select(-rowid)


humanized_data = read_csv("data/developability-humanised.csv") %>% 
  select(params) %>%
  scale_data()

therapeutic_data = read_csv("data/developability-thera.csv") %>% 
  filter(AbOrig_chain == 'heavy' & AbOrig_species == 'human') %>% 
  select(params) %>%
  na.omit() %>% 
  scale_data()

patented_data = read_csv("data/developability-patented.csv") %>% 
  filter(AbOrig_chain == 'heavy' & AbOrig_species == 'human') %>% 
  select(params) %>%
  scale_data()

# import training indices
train_splits = np$load("data/human-heavy_train_splits.npz")

sample_sizes = c('50', '100', '500', '1000', '10000', '20000')


# coefficient of determination formula
rsq <- function (observed, predicted) {
  1 - (sum((observed - predicted)^2) / sum((observed - mean(observed))^2))
}

# Define DP subsets for training
set.seed(69)
dp_sets = list(comb_parameters = params,
               seq_parameters = params[str_detect(params, 'AbChain')],
               struc_parameters = params[str_detect(params, 'AbStruc')],
               comb_parameters_mwds = AbChStruc_mwds,
               seq_parameters_mwds = AbChain_mwds,
               struc_parameters_mwds = AbStruc_mwds,
               comb_parameters_no_mwds = params[!params %in% AbChStruc_mwds],
               seq_parameters_no_mwds = params[str_detect(params, 'AbChain') & !params %in% AbChain_mwds],
               struc_parameters_no_mwds = params[str_detect(params, 'AbStruc') & !params %in% AbStruc_mwds],
               comb_parameters_subsample = sample(params, length(AbChStruc_mwds)),
               seq_parameters_subsample = sample(params[str_detect(params, 'AbChain')], length(AbChain_mwds)),
               struc_parameters_subsample = sample(params[str_detect(params, 'AbStruc')], length(AbStruc_mwds)),
               subsample_mwds = c(sample(AbStruc_mwds, 18),
                                  sample(AbChain_mwds, 10))
)

dp_set_names = names(dp_sets)




# Rest is for reference ######
models = list()
predictions = tibble(parameter = character(),
                     rsq = numeric(),
                     dp_set_name = character(),
                     test_set = character())
sample_size = '1000'
for (dp_set in names(dp_sets)) {
  for (i in 1:20) {
    train_idxs = train_splits$f[[sample_size]][i, ]
    train = AbChStruc_all %>%
      filter(rowid %in% train_idxs) %>%
      select(any_of(dp_sets[[dp_set]]))
    params = train %>% colnames()
    for (param in params) {
      name = glue('{dp_set}_{param}_{sample_size}_{i}')
      
      print(glue('Training {dp_set}_{param}_{sample_size}_{i}'))
      models[[name]] = lm(train[[param]] ~ ., data = train %>% select(-param))

      print(glue('Predicting {dp_set}_{param}_{sample_size}_{i}'))
      predicted_native = predict(models[[name]] , newdata = native_data)
      predicted_humanized = predict(models[[name]] , newdata = humanized_data)
      predicted_therapeutic = predict(models[[name]] , newdata = therapeutic_data)
      predicted_patented = predict(models[[name]] , newdata = patented_data)
      
      predictions = predictions %>%
        add_row(
          parameter = param,
          rsq = rsq(predicted_native, native_data[[param]]),
          dp_set_name = dp_set,
          test_set = 'native'
        ) %>%
        add_row(
          parameter = param,
          rsq = rsq(predicted_humanized, humanized_data[[param]]),
          dp_set_name = dp_set,
          test_set = 'humanized'
        ) %>%
        add_row(
          parameter = param,
          rsq = rsq(predicted_therapeutic, therapeutic_data[[param]]),
          dp_set_name = dp_set,
          test_set = 'therapeutic'
        ) %>%
        add_row(
          parameter = param,
          rsq = rsq(predicted_patented, patented_data[[param]]),
          dp_set_name = dp_set,
          test_set = 'patented'
        )
    }
  }
}
saveRDS(predictions, 'data/predictions_non_native.RDS')
