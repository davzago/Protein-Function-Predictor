from dataset_utils import *
from parsers import *
import argparse
import os
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score, precision_score, recall_score
from skmultilearn.problem_transform import BinaryRelevance
from skmultilearn.model_selection import iterative_train_test_split
import numpy as np
from joblib import dump

parser = argparse.ArgumentParser(description='Training a logistic regressor on a dataset')
parser.add_argument('interpro_set', help='File containing set of possible interpro ids')
parser.add_argument('ontology_file', help='File containing the ontology')
parser.add_argument('dataset_file', help='File containing the dataset for the training of the logistic regression model')
parser.add_argument('-output_path', help='Path to the folder where the predction file will be put', default="InterPro")
args = parser.parse_args()

interpro_set = args.interpro_set
ontology_file = args.ontology_file
dataset_file = args.dataset_file
output_path = args.output_path

if not os.path.isdir(output_path):
    os.mkdir(output_path)

ip_set = parse_interpro_set(interpro_set)
ontology = parse_ontology(ontology_file)
X, y, ip_dict, go_dict, prot_dict = parse_dataset(dataset_file, ip_set, ontology)


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42, shuffle=True)
#X_train, X_test, y_train, y_test = iterative_train_test_split(X, y, test_size=0.33)

y_train, y_test, go_dict = remove_unused_go_terms(y_train, y_test,  go_dict)
save_dict(go_dict, output_path, "go_dict")
save_prot_dict(prot_dict, output_path, "protein_dict")
#save_dict(ip_dict, output_path, "ip_dict")
#print("starting training")
clf = BinaryRelevance(LogisticRegression(random_state=42, verbose=0, solver='saga', max_iter=100)).fit(X_train, y_train)
#print("end training")

dump(clf, output_path + '/InterPro_regression_model.joblib')

y_pred = clf.predict(X_test)
#print("y_pred")
#y_prob = clf.predict_proba(X_test[:20,:])
#print("y_prob")

precision = precision_score(y_test,y_pred, average='micro')
recall = recall_score(y_test,y_pred, average='micro')
f1 = f1_score(y_test,y_pred, average='micro')

with open(output_path + "/metrics.txt", "w") as f:
    f.write("Precision score: " + str(precision) + '\n')
    f.write("Recall score: " + str(recall) + '\n')
    f.write("F1 score: " + str(f1) + '\n')


