from numpy.core.records import array
from parsers import parse_dict, parse_dict_to_array, parse_features
from dataset_utils import *
from dataset_utils import parse_interpro_set, reverse_dict, save_prediction, set_ip_indices
from sklearn.multioutput import MultiOutputClassifier
from sklearn.linear_model import LogisticRegression
import argparse
import numpy as np
import os
from joblib import load


parser = argparse.ArgumentParser(description='Predictor based on InterPro')
parser.add_argument('interpro_set', help='File containing set of possible interpro ids')
parser.add_argument('ontology_dict_file', help='File containing the go terms and their corresponding index') 
parser.add_argument('model', help='Path to the model used to output the predictions')
parser.add_argument('protein_features', help='File containing the proteins to predict and their corresponding features')
parser.add_argument('-output_path', help='Path to the folder where the predction file will be put', default="InterPro_predictions")
args = parser.parse_args()

interpro_set = args.interpro_set
features = args.protein_features
terms_dict = args.ontology_dict_file
model = args.model
output_path = args.output_path

if not os.path.isdir(output_path):
    os.mkdir(output_path)

clf = load(model)
ip_set = parse_interpro_set(interpro_set)
#go_dict = parse_dict(terms_dict)
#go_array = reverse_dict(go_dict)
go_array = parse_dict_to_array(terms_dict)

X, ip_dict, prot_dict = parse_features(features, ip_set)

pred = clf.predict_proba(X)

prot_array = reverse_dict(prot_dict)

Y = reshape_prediction(pred)

pred_dict = prediction_to_dict(Y, prot_array, go_array)

save_prediction(pred_dict, 200, output_path)



