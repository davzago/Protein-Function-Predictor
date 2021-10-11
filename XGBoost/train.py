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
parser.add_argument('goa_file', help="File containing the goa database")
parser.add_argument('k', help="maximum number of proposed go term from each component")
parser.add_argument('-ref_file', help="path to the reference file, if this parameter is set the goa will not be used")# will be substituted with the goa file
parser.add_argument('-model_name', help='Path to the folder where the predction file will be put', default="xgb")
parser.add_argument('-output_path', help='Path to the folder where the predction file will be put', default="output")
args = parser.parse_args()

naive_folder = args.Naive_predictions
blast_folder = args.Blast_predictions
interpro_folder = args.InterPro_predictions
goa_file = args.goa_file
ref = args.ref_file
output_path = args.output_path
k = int(args.k)
model_name = args.model_name

# building the  
naive_pred = parse_component_prediction(naive_folder)
blast_pred = parse_component_prediction(blast_folder)
interpro_pred = parse_component_prediction(interpro_folder)

pred_dict = combine_dictionaries(naive_pred, blast_pred, interpro_pred, k)
#key = random.choice(list(pred_dict))
if ref is not None:
    pred_dict = add_ground_truth(ref, pred_dict)
else:
    pred_dict = add_ground_truth_from_goa(goa_file, pred_dict)

df, assoc = build_dataset(pred_dict)


params = {'objective': 'rank:pairwise', 'learning_rate': 0.01, 'max_depth': 4,
             'gamma': 1.0, 'min_child_weight': 0.1, 'n_estimators': 1000}

X_train = df[['Naive', 'Blast', 'InterPro']].to_numpy().astype('float')
y_train = df['Label'].to_numpy().astype('int')

groups = df.groupby('Group').size().to_frame('size')['size'].to_numpy()

rank_model = xgb.sklearn.XGBRanker(**params)    

rank_model.fit(X_train, y_train, group=groups, verbose=True)

rank_model.save_model(output_path + '/' + model_name + ".model")
