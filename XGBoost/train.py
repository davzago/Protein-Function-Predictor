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
pred_dict = build_prediction_dict([naive_folder,blast_folder,interpro_folder], 30)
#key = random.choice(list(pred_dict))
pred_dict = add_ground_truth(ref, pred_dict)

df, assoc = build_dataset(pred_dict)


params = {'objective': 'rank:pairwise', 'learning_rate': 0.01, 'max_depth': 4, 'gamma': 1.0, 'min_child_weight': 0.1, 'n_estimators': 1000}