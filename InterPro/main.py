from pandas.core.arrays.sparse import dtype
from dataset_utils import *
from parsers import *
from sklearn.multioutput import MultiOutputClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import precision_score
from sklearn.model_selection import train_test_split
from skmultilearn.problem_transform import BinaryRelevance
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix

ip_set = parse_interpro_set("/home/davide/Documenti/Protein Function Predictor/InterPro/interpro_set.txt")
ontology = parse_ontology("/home/davide/Documenti/training/gene_ontology_edit.obo.2018-08-01")
X, y, ip_dict, go_dict, prot_dict = parse_dataset("/home/davide/Documenti/Protein Function Predictor/InterPro/dataset.txt", ip_set, ontology)


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42, shuffle=True)
n = 19742

arr = np.concatenate((np.ones((int(n/2), 5), dtype='bool'),np.zeros((int(n/2),5), dtype='bool')), axis=0)
arr2 = np.concatenate((arr,np.zeros((n,1), dtype=bool)), axis=1)
"""print(y_train)
y_train = csc_matrix(arr2)
print(y_train)

print(X_train.shape, y_train.shape)"""
#data = pd.DataFrame(y_train.toarray())

#print(data)

"""for i in range(0,y_test.shape[1]):
    asd = (y_test.getcol(i))"""

"""row = [1, 2, 3]
col = [0,1,3]
data = [1, 1, 1]
prova = csc_matrix((data, (row,col)), shape=(X_test.shape[0],5))"""

"""matr = np.ones((3,5), dtype='bool')
rest = np.zeros((X_test.shape[0]-3,5), dtype='bool')
matr = np.concatenate((matr, rest), axis=0)
prova2 = csc_matrix(matr)"""

#print(X_train.toarray().any(), y_train[:500,:].toarray().any())

y_train, go_dict = remove_unused_go_terms(y_train[:500,:], go_dict)
print(y_train.shape)
print("starting training")
clf = BinaryRelevance(LogisticRegression(random_state=42)).fit(X_train[:500,:], y_train)
print("end training")

y_pred = clf.predict(X_test)

print(precision_score(y_test, y_pred, average='micro'))





"""prot_index = 0
ip_ids = set()
go_ids = set()

for k, v in prot_dict.items():
    if v == prot_index:
        print(k)
        break

for ip in ip_set:
    ip_idx = ip_dict[ip]
    if x_train[prot_index, ip_idx]:
        ip_ids.add(ip)

for go in ontology:
    go_idx = go_dict[go]
    if y_train[prot_index, go_idx]:
        go_ids.add(go)"""

#print(x_train.shape, y_train.shape)


