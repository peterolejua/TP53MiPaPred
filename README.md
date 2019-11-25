# Tp53MiPaPred

## TP53 missense variants pathogenicity prediction using a Random Forest Meta-predictor

This Github contains the reproducible analysis used for the development of TP53MiPaPred. The following lines describe folders' content and how to predict with this in-silico tool.

### 0_RawData and 1_processed_data

Two datasets were used for the training of Tp53MiPaPred. TP53 missense variants are retrieved from the UMD TP53 database 2017_R2 version. This was merged with the dbNSFP4.0b1a database to retrieve 25 in-silico tools which integrate the new meta-predictor. This folder contains the R project `get_data` and the rmd file `0_get_data2`. These are used to download the raw data mentioned above. Raw data is heavy enough that it was not included in this repo. However, the code fully explain how to download them. The merged data is then saved into the folder 1_processed_data as `data.rds`.

### 1_exploratory_DA

Exploratory data analysis was performed inside this folder and it is available in the html file `1_exploratory_DA`.

### 2_Random_forest

The Random Forest (RF) is trained and tuned inside this folder after imputing data. The rmd file `2-random_forest` develops the whole training. First the problem is introduced. Then missing data is imputed with random forest proximity approach. Parameters for the RF were tuned with a recent methodology implemented with the `tuneRanger` function.

With the tuned parameters an RF was trained using the `ranger` package. After performing a ROC analysis with the OOB data, the RF is modified to find two new models depending on the results found in the pathogenic class. One is meant to be used in settings where more sensitivity is needed (TP53MiPaPred_sens). The second is meant to be used in settings where a higher positive predictive value is needed (TP53MiPaPred_ppv). The new model is applied to Variants of Unknown Significance that were complete cases with the 25 predictors.

#### Results

The imputed data used for training was exported to the file `data_model_imputed.xlsx`. The `model` obtained and the function used to predict with it `TP53MiPaPred(.)` are available as R objects in the file `load_to_predict.RData`. `model$predictions` are the OOB group membership probabilities. These can be used to further tuned TP53MiPaPred in case sensibility or PPV for the pathogenic class had to be modified. Finally, the user will require a `newdata` with complete cases in the 25 predictors to predict with this model.

## 3_thesis_document

This folder contains a pdf file that explains a thesis document corresponding to the whole development of TP53MiPaPred.

