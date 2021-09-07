from utils import *
import argparse
import pandas as pd
import xgboost as xgb
from sklearn.utils import shuffle
import sklearn as sk

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
pred_dict = build_prediction_dict([naive_folder,blast_folder,interpro_folder]) 
pred_dict = add_ground_truth(ref, pred_dict)

df, assoc = build_dataset(pred_dict)

X = df[['Naive', 'Blast', 'InterPro']].to_numpy()
y = df['Label'].to_numpy()

print(X)

n_rows = len(df.index)
split = int(n_rows * 0.8)

#train_rows = df.iloc[:, :split]
#test_rows = df.iloc[:, split:]

params = {'objective': 'rank:pairwise', 'learning_rate': 0.1,
          'gamma': 1.0, 'max_depth': 6}

"""gss = GroupShuffleSplit(test_size=.40, n_splits=5, random_state = 7).split(df, groups=df['id'])

X_train_inds, X_test_inds = next(gss)

train_data= df.iloc[X_train_inds]
X_train = train_data.loc[:, ~train_data.columns.isin(['id','rank'])]
y_train = train_data.loc[:, train_data.columns.isin(['rank'])]

groups = train_data.groupby('id').size().to_frame('size')['size'].to_numpy()

test_data= df.iloc[X_test_inds]

#We need to keep the id for later predictions
X_test = test_data.loc[:, ~test_data.columns.isin(['rank'])]
y_test = test_data.loc[:, test_data.columns.isin(['rank'])]"""

kFold = GroupKFold(n_splits=5)

scores = []


for train_index, test_index in kFold.split(df[['Naive', 'Blast', 'InterPro']], df['Labels'], groups=df['Groups']):

    train_data = df.iloc[train_index]

    X_train = train_data[['Naive', 'Blast', 'InterPro']]
    y_train = train_data['Labels']

    groups = train_data.groupby('Group').size().to_frame('size')['size'].to_numpy()

    test_data = df.iloc[test_index]

    X_test = test_data[['Naive', 'Blast', 'InterPro']]
    y_test = test_data['Label']
    
    rank_model = xgb.sklearn.XGBRanker(**params)    

    rank_model.fit(X_train, y_train, verbose=True)
    pred = rank_model.predict(X_test)

    f1_score = sk.f1_score(y_test, pred)
    scores.append(f1_score)

print(scores)



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






