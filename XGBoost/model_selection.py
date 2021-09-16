from utils import *
import argparse
import pandas as pd
import xgboost as xgb
from sklearn.utils import shuffle
import sklearn as sk
import random

from sklearn.model_selection import cross_val_score, GridSearchCV, GroupKFold, GroupShuffleSplit


parser = argparse.ArgumentParser(description='XGBoost predictor')
parser.add_argument('Naive_predictions', help='Path to the folder containing the predictions from the Naive component')
parser.add_argument('Blast_predictions', help='Path to the folder containing the predictions from the Blast component')
parser.add_argument('InterPro_predictions', help="Path to the folder containing the predictions from the Interpro component")
parser.add_argument('ref_file', help="path to the reference file")# will be substituted with the goa file
parser.add_argument('-output_path', help='Path to the folder where the predction file will be put', default="output")
args = parser.parse_args()

naive_folder = args.Naive_predictions
blast_folder = args.Blast_predictions
interpro_folder = args.InterPro_predictions
ref = args.ref_file

naive_pred = parse_component_prediction(naive_folder)
blast_pred = parse_component_prediction(blast_folder)
interpro_pred = parse_component_prediction(interpro_folder)

pred_dict = combine_dictionaries(naive_pred, blast_pred, interpro_pred, 10)

#pred_dict = build_prediction_dict([naive_folder,blast_folder,interpro_folder], 50)
key = random.choice(list(pred_dict))
print(key)
print(pred_dict[key])

pred_dict = add_ground_truth(ref, pred_dict)

df, assoc = build_dataset(pred_dict)


params = {'objective': 'rank:pairwise', 'learning_rate': 0.01, 'max_depth': 4, 'gamma': 1.0, 'min_child_weight': 0.1, 'n_estimators': 1000}
"""params = {'objective': 'rank:pairwise', 'learning_rate': 0.1,
          'gamma': 1.0, 'min_child_weight': 0.1,
          'max_depth': 6, 'n_estimators': 4}"""

gs_params = {'objective': ['rank:pairwise'], 'gamma': [1.0, 10.0, 0.1, 0.01], 
                'max_depth': [6,4,10], 'min_child_weight': [0.1],
                'n_estimators': [10, 100, 1000]}# 'tau': [0.3, 0.4, 0.6, 0.8]}

objective, gamma, max_depth, min_child_weight, n_estimators = list(gs_params.values())
params = [[obj, g, depth, weight, estim] for obj in objective for g in gamma for depth in max_depth for weight in min_child_weight for estim in n_estimators]
       
kFold = GroupKFold(n_splits=5)

mean = []
parameters_list = []
for obj, g, d, w, e  in params:
    f1_scores = []
    precision_scores = []
    recall_scores = []
    parameters = {"objective": obj, "gamma": g, 'max_depth': d, "min_child_weight": w, "n_estimators": e}
    parameters_list.append(parameters)

    for train_index, test_index in kFold.split(df[['Naive', 'Blast', 'InterPro']], df['Label'], groups=df['Group']):
        
        train_data = df.iloc[train_index]

        X_train = train_data[['Naive', 'Blast', 'InterPro']].to_numpy().astype('float')
        y_train = train_data['Label'].to_numpy().astype('int')

        groups = train_data.groupby('Protein id').size().to_frame('size')['size'].to_numpy()

        test_data = df.iloc[test_index]

        X_test = test_data[['Naive', 'Blast', 'InterPro']].to_numpy().astype('float')
        y_test = test_data['Label'].to_numpy().astype('bool')

        test_groups = test_data.groupby('Group').size().to_frame('size')['size'].to_numpy()
        
        rank_model = xgb.sklearn.XGBRanker(**parameters)    

        rank_model.fit(X_train, y_train, group=groups, verbose=True)

        #pred = X_test.groupby('Group').apply(lambda x: predict(rank_model, x))
        """pred = rank_model.predict(X_test)
        print((pred > 1).any())"""
        #print(y_train)
        pred = predict2(rank_model, X_test, test_groups)

        tau = 0.8

        f1_score = sk.metrics.f1_score(y_test, (pred > tau))
        precision_score = sk.metrics.precision_score(y_test, (pred > tau))
        recall_score = sk.metrics.recall_score(y_test, (pred > tau))

        f1_scores.append(f1_score)
        precision_scores.append(precision_score)
        recall_scores.append(recall_score)
    
    print("Parameters:", parameters)
    print("F1:", f1_scores)
    print("Precision:", precision_scores)
    print("Recall:", recall_scores)
    mean.append(scores_mean(f1_scores))
    print()

print("the means are:")
print(mean)
print("the best parameters based on F score are:")
max_val = max(mean)
max_idx = mean.index(max_val)
best_params = parameters_list[max_idx]
print(best_params)
print("with mean:", max_val)


# tau = 0 ; k = 10 ; non normalized scores
# the best parameters based on F score are:
# {'objective': 'rank:pairwise', 'gamma': 10.0, 'max_depth': 10, 'min_child_weight': 0.1, 'n_estimators': 1000}
# with mean: 0.4741091645761764

# tau = 0.5 ; k = 10 ; non normalized scores
# the best parameters based on F score are:
# {'objective': 'rank:pairwise', 'gamma': 0.1, 'max_depth': 4, 'min_child_weight': 0.1, 'n_estimators': 100}
# with mean: 0.5001041800573793

# tau = 0.5 ; k = 30 ; non normalized scores
# the best parameters based on F score are:
# {'objective': 'rank:pairwise', 'gamma': 0.1, 'max_depth': 10, 'min_child_weight': 0.1, 'n_estimators': 100}
# with mean: 0.38482721885209653

# tau = 0.7 ; k = 30 ; normalized scores
# the best parameters based on F score are:
# {'objective': 'rank:pairwise', 'gamma': 10.0, 'max_depth': 4, 'min_child_weight': 0.1, 'n_estimators': 100}
# with mean: 0.38793471738974544

# tau = 0.7 ; k = 10 ; normalized scores
# the best parameters based on F score are:
# {'objective': 'rank:pairwise', 'gamma': 10.0, 'max_depth': 4, 'min_child_weight': 0.1, 'n_estimators': 10}
# with mean: 0.4994602001762584

# tau = 0.5 ; k = 50 ; normalized scores
# the best parameters based on F score are:
# {'objective': 'rank:pairwise', 'gamma': 1.0, 'max_depth': 4, 'min_child_weight': 0.1, 'n_estimators': 10}
# with mean: 0.3070115179051105

# tau = 0.8 ; k = 50 ; normalized scores
# the best parameters based on F score are:
# {'objective': 'rank:pairwise', 'gamma': 10.0, 'max_depth': 4, 'min_child_weight': 0.1, 'n_estimators': 10}
# with mean: 0.3448914151305004


