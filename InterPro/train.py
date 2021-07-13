from dataset_utils import *
from parsers import *
import argparse
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import precision_score
from sklearn.model_selection import train_test_split
from skmultilearn.problem_transform import BinaryRelevance
from sklearn.multioutput import MultiOutputClassifier
import numpy as np
from joblib import dump

parser = argparse.ArgumentParser(description='dataset builder')
parser.add_argument('interpro_set', help='File containing set of possible interpro ids')
parser.add_argument('ontology_file', help='File containing the ontology')
parser.add_argument('dataset_file', help='File containing the dataset for the training of the logistic regression model')
parser.add_argument('-output_path', help='Path to the folder where the predction file will be put', default="output")
args = parser.parse_args()

interpro_set = args.interpro_set
ontology_file = args.ontology_file
dataset_file = args.dataset_file
output_path = args.output_path

ip_set = parse_interpro_set(interpro_set)
ontology = parse_ontology(ontology_file)
X, y, ip_dict, go_dict, prot_dict = parse_dataset(dataset_file, ip_set, ontology)


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42, shuffle=True)

y_train, y_test, go_dict = remove_unused_go_terms(y_train[:20,:], y_test,  go_dict)
print(y_train.shape)
print("starting training")
clf = BinaryRelevance(LogisticRegression(random_state=42, verbose=0, solver='saga', max_iter=30)).fit(X_train[:20,:], y_train)
print("end training")

dump(clf, 'InterPro_regression_model_2.joblib')

"""clf = load('InterPro_regression_model.joblib')
"""
y_pred = clf.predict(X_test[:100,:])

probs = clf.predict_proba(X_test[:100,:])

print(probs)


print(precision_score(y_test[:100], y_pred, average='micro'))
