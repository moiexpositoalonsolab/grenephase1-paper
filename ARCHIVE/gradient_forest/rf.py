import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
import pickle
import dask.dataframe as dd
import seaborn as sns
import numpy as np

from sklearn.metrics import r2_score

with open('../leave_1_out/splits_samples.pkl', 'rb') as file:
    splits_samples = pickle.load(file)

train_samples = splits_samples[0][0]
test_samples = splits_samples[0][1]

## import deltap ld for testing because not ldp will take forever 
delta_ldp = pd.read_csv('/carnegie/nobackup/scratch/xwu/grenet/merged_frequency/merged_hapFIRE_delta_p_LDpruned.txt', 
                        sep = '\t', 
                        usecols = train_samples)

clim_sites_during_exp = pd.read_csv('/carnegie/nobackup/scratch/tbellagio/grene/data/bioclimvars_experimental_sites_era5.csv')

sites = pd.Series(train_samples).str.split('_').str[0].astype(int)

env = sites.reset_index().merge(clim_sites_during_exp, right_on = 'site', left_on = 0).drop(['index'],axis=1)
env = env.drop(0,axis=1)
mean_env = env['bio1'].mean()
std_env = env['bio1'].std()

env['bio1'] = (env['bio1'] - mean_env) / std_env

# Load data
#data = pd.read_csv('data.csv')

# Define the features and targets
X_train = np.array(env.drop('site',axis=1).copy()['bio1']).reshape(-1, 1)#.reset_index()  # n_env_vars is the number of environmental variables
#X_train = env.drop('site',axis=1).copy()
y_train_full = delta_ldp.T.copy()

sites_test = pd.Series(test_samples).str.split('_').str[0].astype(int)

env_test = sites_test.reset_index().merge(clim_sites_during_exp, right_on = 'site', left_on = 0).drop(['index'],axis=1)
env_test = env_test.drop(0,axis=1)
#X_test = env_test.drop('site',axis=1).copy()

env_test['bio1'] = (env_test['bio1'] - mean_env) / std_env

X_test = np.array(env_test.drop('site',axis=1).copy()['bio1']).reshape(-1, 1)

## import deltap ld for testing because not ldp will take forever 
delta_ldp_test = pd.read_csv('/carnegie/nobackup/scratch/xwu/grenet/merged_frequency/merged_hapFIRE_delta_p_LDpruned.txt', 
                        sep = '\t', 
                        usecols = test_samples)

y_test_full = delta_ldp_test.T.copy()



from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import make_scorer, r2_score

param_grid = {
    'n_estimators': [100, 200, 300],  # Number of trees in the forest
    #'max_features': ['auto', 'sqrt'],  # Number of features to consider at every split
    'max_depth': [2, 3, 5, 10, 20],   # Maximum number of levels in tree
     # Minimum number of samples required at each leaf node
}

# Create a random forest regressor object
forest = RandomForestRegressor(random_state=42)

# Setup the grid search with cross-validation
grid_search = GridSearchCV(estimator=forest, param_grid=param_grid, 
                           scoring='r2', cv=3, verbose=2, n_jobs=-1)

r2_scores_train_all = []
import pandas as pd

r2_scores_train_all = []
r2_scores_test_all = []
best_params_all = []

import os

# Open a text file for appending results. This will create the file if it doesn't exist.
if not os.path.exists('model_results.txt'):
    with open('model_results.txt', 'w') as file:
        file.write("SNP,Best Parameters,R2 Score Train,R2 Score Test\n")

for i in y_train_full.columns:
    y_train = y_train_full[i]
    y_test = y_test_full[i]
    
    # Fit the grid search
    grid_search.fit(X_train, y_train)
    
    # Best estimator and parameters
    best_forest = grid_search.best_estimator_
    best_params = grid_search.best_params_
    
    # Evaluate on the training set
    y_predict_train = best_forest.predict(X_train)
    
    X = X_train.reshape(-1, 1)  # Independent variable
    y = y_predict_train # Dependent variable

    # Fit the linear model
    model = LinearRegression().fit(X, y)
    # Predict the values
    y_pred = model.predict(X)
    
    r2_score_train = r2_score(y, y_pred)

    #r2_score_train = r2_score(y_train, y_predict_train)
    
    # Evaluate on the testing set
    y_predict_test = best_forest.predict(X_test)


    X = X_test.reshape(-1, 1)  # Independent variable
    y = y_predict_test # Dependent variable

    # Fit the linear model
    model = LinearRegression().fit(X, y)
    # Predict the values
    y_pred = model.predict(X)
    
    r2_score_test = r2_score(y, y_pred)
    
    #r2_score_test = r2_score(y_test, y_predict_test)
    
    # Format the best parameters for easier reading in the text file
    params_formatted = ', '.join(f"{k}: {v}" for k, v in best_params.items())
    
    # Append results to the text file
    with open('model_results.txt', 'a') as file:
        file.write(f"{i},{params_formatted},{r2_score_train},{r2_score_test}\n")
    
    # Print best parameters and R2 score for this SNP
    print(f"Best parameters for SNP {i}: {best_params}")
    print(f"R2 score on test set for SNP {i}: {r2_score_test}")

