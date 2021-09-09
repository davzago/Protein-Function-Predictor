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

# building the  
pred_dict = build_prediction_dict([naive_folder,blast_folder,interpro_folder],70)
print(len(random.choice(list(pred_dict.values()))))
key = random.choice(list(pred_dict))
pred_dict = add_ground_truth(ref, pred_dict)

df, assoc = build_dataset(pred_dict)


params = {'objective': 'rank:pairwise', 'learning_rate': 0.01, 'max_depth': 4, 'gamma': 1.0, 'min_child_weight': 0.1, 'n_estimators': 1000}
"""params = {'objective': 'rank:pairwise', 'learning_rate': 0.1,
          'gamma': 1.0, 'min_child_weight': 0.1,
          'max_depth': 6, 'n_estimators': 4}"""

gs_params = {'objective': ['rank:pairwise'], 'gamma': [1.0, 10.0, 0.1, 0.01], 
                'max_depth': [6,4,10], 'min_child_weight': [0.1],
                'n_estimators': [10, 100, 1000]}

kFold = GroupKFold(n_splits=5)

f1_scores = []
precision_scores = []
recall_scores = []

"""X = df[['Naive', 'Blast', 'InterPro']].to_numpy().astype('float')
y = df['Label'].to_numpy().astype('int')
groups = df.groupby('Protein id').size().to_frame('size')['size'].to_numpy()

rank_model = xgb.sklearn.XGBRanker()

ranker = GridSearchCV(estimator=rank_model, param_grid=gs_params)

ranker.fit(X, y, groups=groups)"""


for train_index, test_index in kFold.split(df[['Naive', 'Blast', 'InterPro']], df['Label'], groups=df['Group']):

    train_data = df.iloc[train_index]

    X_train = train_data[['Naive', 'Blast', 'InterPro']].to_numpy().astype('float')
    y_train = train_data['Label'].to_numpy().astype('int')

    groups = train_data.groupby('Protein id').size().to_frame('size')['size'].to_numpy()

    test_data = df.iloc[test_index]

    X_test = test_data[['Naive', 'Blast', 'InterPro']].to_numpy().astype('float')
    y_test = test_data['Label'].to_numpy().astype('bool')

    test_groups = test_data.groupby('Group').size().to_frame('size')['size'].to_numpy()
    
    rank_model = xgb.sklearn.XGBRanker(**params)    

    rank_model.fit(X_train, y_train, group=groups, verbose=True)

    #pred = X_test.groupby('Group').apply(lambda x: predict(rank_model, x))
    """pred = rank_model.predict(X_test)
    print((pred > 1).any())"""
    #print(y_train)
    pred = predict2(rank_model, X_test, test_groups)
    print((pred > 1).any())
    print(pred)

    tau = 0.8

    f1_score = sk.metrics.f1_score(y_test, (pred > tau))
    precision_score = sk.metrics.precision_score(y_test, (pred > tau))
    recall_score = sk.metrics.recall_score(y_test, (pred > tau))

    f1_scores.append(f1_score)
    precision_scores.append(precision_score)
    recall_scores.append(recall_score)


print("F1:", f1_scores)
print("Precision:", precision_scores)
print("Recall:", recall_scores)



#train = xgb.DMatrix(train_rows['Naive', 'Blast', 'InterPro'], label=train_rows['Label'])
#train.set_group(train_rows['group'])

#test = xgb.DMatrix(test_rows['Naive', 'Blast', 'InterPro'], label=test_rows['Label'])
#test.set_group(test_rows['group'])

#x_train = train_rows['Naive', 'Blast', 'InterPro']
#y_train = train_rows['Label']
#group_train = train_rows['Group']

#x_test = test_rows['Naive', 'Blast', 'InterPro']
#y_test = test_rows['Label']
#group_test = test_rows['Group']






