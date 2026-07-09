library(tidyverse)
library(ranger)
library(caret)
library(doFuture)
library(futurize)
source('/Gages2 Data/RF model/utilities.R')

#settings for parallelization in the cluster
# slurm_cores = as.numeric(Sys.getenv('cores'))
# print(paste('CPUs assigned by SLURM:',slurm_cores))
# plan(multicore, workers = slurm_cores)

#choose which signatures to model
metrics_sel <- c('Q_mean', 'Q95', 'Q10', 'Q95_minus_Q10', 'TotalRR', 'BFI',
                 'qp_elasticity', 'Q_frequency_high_3', 'Q_frequency_noflow',
                 'Q_totalduration_high_3', 'Q_meanduration_high_3',
                 'Q_totalduration_noflow', 'HFD_mean', 'HFI_mean',
                 'peakQ_timing', 'FlashinessIndex', 'RLD', 'FDC_slope',
                 'BaseflowRecessionK', 'Recession_a_Seasonality')

#load in all the data
rf_data <- read_csv('data/gages/predictors/rf_model_data.csv')

#select which predictors to use
catchment_predictors = c('precip_annual', 'pet_annual',
                         'precip_jfm', 'pet_jfm',
                         'precip_amj', 'pet_amj',
                         'precip_jas', 'pet_jas',
                         'precip_ond', 'pet_ond',
                         'si', 'swe_annual', 'max_swe',
                         'zero_swe_day', 'swe_persistence', 'melt_duration',
                         'ag', 'developed', 'forest', 'grass',
                         'elev', 'slope', 'twi',
                         'dist_index', 'water_use_mean',
                         'tile_pct', 'dam_storage',
                         'soil_perm', 'soil_depth', 'age')
n_preds <- length(c(catchment_predictors))

#set up hyperparameter grid for tuning
#number of parameter sets to test. 60 is a good rule of thumb but can bump up to 100 if training is fast
param_searches <- 100
#create grid using set ranges of these three parameters. Can adjust as needed, good reference at:
#https://bradleyboehmke.github.io/HOML/random-forest.html#hyperparameters
param_grid <- expand.grid(mtry = 2:20,
                          min.node.size = 3:30,
                          splitrule = c('variance', 'extratrees'))
#randomly select parameter sets from the grid
set.seed(527)
sample_indices <- sample(seq(1, nrow(param_grid)), size = param_searches)
param_grid_sample <- param_grid[sample_indices,]

#Initialize ranking table
all_ds_ranks <- tibble()

#loop through each metric to train models
for(m in metrics_sel){
  print(paste('Starting model for metric:',m))
  
  #set name for model and filepaths to results
  model_name <- m
  # model_name <- paste0(tolower(m),'_cv5')
  if(!dir.exists(paste0('data/models/',model_name))) dir.create(paste0('data/models/',model_name))
  if(!dir.exists(paste0('data/models/',model_name,'/training'))) dir.create(paste0('data/models/',model_name,'/training'))
  
  #set up all the data for modeling
  #select the metric and all predictors to be used
  response <- paste0('hw_', m)
  ds_metric <- paste0('ds_', m)
  predictors <- c(
    ds_metric,
    paste0('hw_', catchment_predictors),
    paste0('ds_', catchment_predictors),
    'drainage_ratio'  #or ds_drainage_area
  )
  
  dat <- rf_data %>%
    select(
      headwater_id, downstream_id, hw_region, all_of(response), all_of(predictors)
    ) %>%
    mutate(across(everything(),
                  ~ifelse(is.infinite(.x) | is.nan(.x), NA, .x))) %>%
    filter(if_all(everything(), ~!is.na(.x))) %>%
    rename(obs = all_of(response))
  
  #Row indexing for joining predictions later
  dat$rowIndex <- seq_len(nrow(dat))
  
  #Model data with no identifiers
  model_dat <- dat %>%
    select(-headwater_id, -downstream_id, -hw_region, -rowIndex)
  
  #Repeated cross-validation setup
  cv_control <- trainControl(
    method = 'repeatedcv',
    number = 10,
    repeats = 5,
    savePredictions = 'final',
    summaryFunction = defaultSummary,
    verboseIter = T
  )
  
  set.seed(527)
  
  #train models using tuning settings. This function will loop through all of our
  #parameter sets and save the model with the parameters that give the best performance
  rf_fit <- train(obs ~ ., 
                  data = model_dat,
                  method = 'ranger',
                  num.trees = 1000,
                  importance = 'permutation',
                  verbose = T,
                  metric = 'MAE',
                  maximize= F,
                  trControl = cv_control,
                  tuneGrid = param_grid_sample) #|>
    #this enables parallelization when on the cluster
    # futurize()
  
  #save the results of hyperparameter tuning, arranged in descending performance order
  #Top row will give you the optimal parameters to use in final model
  tune_results <- rf_fit$results %>%
    arrange(MAE)
  write_csv(tune_results, paste0('data/models/',model_name,'/training/tune_results.csv'))
  
  #run test data through the best model to assess performance on new data
  test_preds <- rf_fit$pred %>%
    inner_join(rf_fit$bestTune,
              by = c('mtry', 'min.node.size', 'splitrule')) %>%
    left_join(
      dat %>% select(rowIndex, hw_region, headwater_id, downstream_id),
      by = 'rowIndex'
    )
  write_csv(test_preds, paste0('data/models/',model_name,'/training/test_preds.csv'))
  
  #Region-wise performance
  test_perf <- test_preds %>%
    group_by(hw_region) %>%
    summarize(
      MAE = mae(pred, obs),
      NSE = nse(pred, obs),
      n=n(),
      .groups = 'drop'
    ) %>%
    mutate(metric = m)
  write_csv(test_perf, paste0('data/models/',model_name,'/training/test_perf.csv'))
  
  #Variable importance
  var_imp <- varImp(rf_fit)
  importance_df <- var_imp$importance %>%
    rownames_to_column('predictor') %>%
    rename(importance = Overall) %>%
    arrange(desc(importance)) %>%
    mutate(rank = row_number(), importance_pct = 100 * importance / sum(importance))
  write_csv(importance_df, paste0('data/models/',model_name,'/training/variable_importance.csv'))
  
  #Save top 10 predictors
  top10 <- importance_df %>%
    slice_head(n = 10)
  write_csv(top10, paste0('data/models/',model_name,'/training/top10_predictors.csv'))
  
  #Downstream predictor rank
  ds_name <- paste0('ds_', m)
  ds_rank <- importance_df %>%
    filter(predictor == ds_name) %>%
    mutate(metric = m, n_predictors = nrow(importance_df))%>%
    select(metric, predictor, rank, importance, importance_pct, n_predictors)
  all_ds_ranks <- bind_rows(all_ds_ranks, ds_rank)
}

write_csv(all_ds_ranks, 'data/models/downstream_predictor_summary.csv')
