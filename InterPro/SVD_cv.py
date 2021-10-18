import os 
import argparse
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import PrecisionRecallDisplay
from sklearn.decomposition import TruncatedSVD
from sklearn.model_selection import KFold
from sklearn.linear_model import LogisticRegression
from sklearn.multioutput import MultiOutputClassifier
from dataset_utils import *
from parsers import *


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
    
n_proteins = 1000

ip_dict = parse_interpro_list(interpro_set)
ontology = parse_ontology(ontology_file)
X, y, ip_dict, go_dict, prot_dict = parse_dataset(dataset_file, ip_dict, ontology)
#X = X[:n_proteins]
#y = y[:n_proteins]

#y, y_XD, go_dict = remove_unused_go_terms(y, y[:1], go_dict)


kf = KFold(n_splits=5, random_state=42, shuffle=True)
fold = 1
for train_index, val_index in kf.split(X):
    X_train, X_val = X[train_index], X[val_index]
    y_train, y_val = y[train_index], y[val_index]

    y_train, y_val, _ = remove_unused_go_terms(y_train, y_val, go_dict)
    
    lsa = TruncatedSVD(n_components=10000, n_iter=5, random_state=42)
    X_train = lsa.fit_transform(X_train)
    print(X_train.shape)
    #print(lsa.explained_variance_ratio_)
    X_val = lsa.fit_transform(X_val)
    print(X_val.shape)
    #print(lsa.explained_variance_ratio_)


    clf = MultiOutputClassifier(LogisticRegression(random_state=42, verbose=0, solver='saga', max_iter=100), n_jobs=8).fit(X_train, y_train.toarray())
    y_score = clf.predict_proba(X_val)

    y_score = reshape_prediction(y_score)

    precision, recall, _ = precision_recall_curve(y_val.toarray().flatten(), y_score.flatten())
    display = PrecisionRecallDisplay(precision, recall)
    display.plot()
    plt.savefig(output_path + '/' + "fold" + str(fold) + ".png")
    fmax = max((2 * precision * recall) / (precision + recall))
    with open(output_path + '/' + "Fmax.txt", 'a') as f:
        f.write("the Fmax for the fold " + str(fold) + " is: " + str(fmax) + '\n')
    fold += 1


