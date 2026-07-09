# Gages2-RandomForest
This code is for a Random Forest model that will predict 20 hydrologic signatures for headwater-downstream pairs across 7 ecoregions (Appalachian Mountains, Central Plains, Eastern Forests, Great Plains, Mountain West, Northern Forests, Desert Southwest).

To run the code, you need the csv file called rf_model_data.csv which contains all the data (sites, predictors, response variables) which is detailed below.

The paired headwater-downstream sites used can be found in the csv hw_ds_connections_exp.csv.

Predictor variables include climate, land cover, storage/soils, topography, water use, and other.  Below is a list of the predictors used in each of the 6 categories.  Each predictor (except the signature metric) is included twice - once for the headwater catchment and once for the downstream catchment:
- Climate - precip_annual, pet_annual, precip_jfm → precip in january, february, march, pet_jfm → pet in january, february, march, precip_amj → precip in april, may, june, pet_amj → pet in april, may, june, precip_jas → precip in july, august, september, pet_jas → pet in july, august, september, precip_ond → precip in october, november, december, pet_ond → pet in october, november, december, swe_annual, max_swe, zero_swe_day, swe_persistence, melt_duration
- Land Cover - ag, developed, forest, grass
- Storage/Soils - si, dist_index, dam_storage, soil_perm, soil_depth
- Topography - elev, slope, twi
- Water Use - water_use_mean
- Other - tile_pct, age, Corresponding downstream signature value


Additionally, the RF model is run using either drainage_ratio which is the drainage area ratio between the paired headwater-downstream catchments or ds_drainage_ratio which is the paired downstream drainage area.

For every hydrologic signature, the response variable is the observed headwater signature value which was predicted using all of the predictors above.

Any missing, undefined (NaN) or infinite values were removed before the model was trained.

Random Forest regression was implemented using Ranger which is a part of the Caret package in R.

Hyperparameter tuning → randomized search of 100 parameter combinations sampled from a grid consisting of:
- Number of variables considered at each split (mtry): 2 - 20
- Minimum terminal node size (min.node.size): 3 - 20
- Split rule: variance or extremely randomized trees (extratrees)

Model performance → used repeated 10-fold cross-validation with 5 repetitions which results in 50 validation folds for each candidate parameter set.

Each RF was trained with 1000 trees.

Mean Absolute Error (MAE) was used as the optimization criterion during the hyperparameter tuning process.  The parameter combination with the lowest MAE was chosen as the final model.

Cross-validation predictions from the chosen final model were stored and used for performance assessments.

Model performance was computed separately for each region using the following:
- Mean Absolute Error (MAE)
- Nash-Sutcliffe Efficiency (NSE)

Predictor importance was calculated using permutation-based variable importance.  Predictor importance scores were ranked for each signature and the relative importance percentages were calculated by normalizing the importance scores to sum to 100%.  The 10 most important predictors for each signature were identified and recorded for further analysis.
