import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
import sys
import pickle
import dask.dataframe as dd
import seaborn as sns
import numpy as np
import os
from sklearn.metrics import r2_score

# Parse command line arguments
number = int(sys.argv[1])
train_samples = sys.argv[2].strip('"').split()
test_samples = sys.argv[3].strip('"').split()


# File and directory setup
directory_name = f"split_number_{number}"
os.makedirs(directory_name, exist_ok=True)  # Creates the directory if it doesn't exist



## import deltap ld for testing because not ldp will take forever 
delta_ldp = pd.read_csv('/carnegie/nobackup/scratch/xwu/grenet/merged_frequency/merged_hapFIRE_delta_p_LDpruned.txt', sep = '\t')

clim_sites_during_exp = pd.read_csv('/carnegie/nobackup/scratch/tbellagio/grene/data/bioclimvars_experimental_sites_era5.csv')

forest = RandomForestRegressor(random_state=42, n_estimators = 100, max_depth = 2)
    
## traning climate 
sites_t = pd.Series(train_samples).str.split('_').str[0].astype(int)
env_t = sites_t.reset_index().merge(clim_sites_during_exp, right_on = 'site', left_on = 0).drop(['index'],axis=1)
env_t = env_t.drop(0,axis=1)
mean_env_t = env_t['bio1'].mean()
std_env_t = env_t['bio1'].std()

env_t['bio1'] = (env_t['bio1'] - mean_env_t) / std_env_t

## testing sites
sites_test = pd.Series(test_samples).str.split('_').str[0].astype(int)
env_test = sites_test.reset_index().merge(clim_sites_during_exp, right_on = 'site', left_on = 0).drop(['index'],axis=1)
env_test = env_test.drop(0,axis=1)
env_test['bio1'] = (env_test['bio1'] - mean_env_t) / std_env_t

# Define the features and targets
X_train = np.array(env_t.drop('site',axis=1).copy()['bio1']).reshape(-1, 1)#.reset_index()  # n_env_vars is the number of environmental variables
#X_train = env.drop('site',axis=1).copy()
y_train_full = delta_ldp[train_samples].T

X_test = np.array(env_test.drop('site',axis=1).copy()['bio1']).reshape(-1, 1)

y_test_full = delta_ldp[test_samples].T

r2_scores_all = {}
predictions = {}
for snp in y_train_full.columns:
    y_train = y_train_full[snp]
    # Train the model
    forest.fit(X_train, y_train)
    ## evalaute on the training 
    y_predict_train = forest.predict(X_train)
    r2_scores_train = r2_score(y_train, y_predict_train)

    if r2_scores_train > 0.2:
        ##evaluate on the testing 
        y_test = y_test_full[snp]
        y_predict_test = forest.predict(X_test)
        r2_scores_test = r2_score(y_test, y_predict_test)
        r2_scores_all[snp] = [r2_scores_train, r2_scores_test]
        predictions[snp] = y_predict_test

r2_scores_all = pd.DataFrame(r2_scores_all).T.reset_index()
r2_scores_all.columns = ['snp_index', 'r2train', 'r2test']
r2_scores_all.to_csv(f'r2/rf_r2_split{number}.csv')

predictions = pd.DataFrame(predictions).T
predictions.columns = test_samples
predictions.to_csv(f'predictions/af_predictions_split{number}.csv')
print(f'done{number}')