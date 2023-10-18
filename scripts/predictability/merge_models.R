library(glue)
library(dplyr)
models = list()
for (i in 1:20) {
  models = models %>% 
    append(readRDS(glue('/storage/jahnz/R/developability/OLS/matteo/model_partitions/parallel_models_{i}.RDS')))
  print(glue('/storage/jahnz/R/developability/OLS/matteo/model_partitions/parallel_models_{i}.RDS'))
}

models = models %>% 
  append(readRDS('/storage/jahnz/R/developability/OLS/matteo/models_1000.RDS'))

saveRDS(models,
        '/storage/jahnz/R/developability/OLS/matteo/all_models.RDS')