from dataset_utils import *
from parsers import *
import argparse
import os
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.multioutput import MultiOutputClassifier
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
parser.add_argument('-model_name', help='name of the output model', default='InterPro_regression_model')
args = parser.parse_args()

interpro_set = args.interpro_set
ontology_file = args.ontology_file
dataset_file = args.dataset_file
output_path = args.output_path
model_name = args.model_name

if not os.path.isdir(output_path):
    os.mkdir(output_path)

ip_dict = parse_interpro_list(interpro_set)
ontology = parse_ontology(ontology_file)
X, y, ip_dict, go_dict, prot_dict = parse_dataset(dataset_file, ip_dict, ontology)


#X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42, shuffle=True)
#X_train, X_test, y_train, y_test = iterative_train_test_split(X, y, test_size=0.33)

y, y_test, go_dict = remove_unused_go_terms(y, y[1,:],  go_dict) # y[1,:] is just to put an argument there, the corresponding return will not be used
#y_train , y_test, go_dict = remove_unused_go_terms(y_train, y_test,  go_dict)
print("Numero di predittori:", len(go_dict))
save_dict(go_dict, output_path, "go_dict")
save_prot_dict(prot_dict, output_path, "protein_dict")
save_dict(ip_dict, output_path, "ip_dict")
#print("starting training")
clf = MultiOutputClassifier(LogisticRegression(random_state=42, verbose=0, solver='saga', max_iter=1000), n_jobs=8).fit(X, y.toarray())
#print("end training")

dump(clf, output_path + '/' + model_name + '.joblib')

#y_pred = clf.predict(X_test)
#print("y_pred")
#y_prob = clf.predict_proba(X_test[:20,:])
#print("y_prob")

"""precision = precision_score(y_test,y_pred, average='micro')
recall = recall_score(y_test,y_pred, average='micro')
f1 = f1_score(y_test,y_pred, average='micro')

with open(output_path + "/metrics.txt", "w") as f:
    f.write("Precision score: " + str(precision) + '\n')
    f.write("Recall score: " + str(recall) + '\n')
    f.write("F1 score: " + str(f1) + '\n')"""



