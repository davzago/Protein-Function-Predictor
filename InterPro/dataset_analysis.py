import os
import argparse
from pandas import DataFrame
import scipy
from dataset_utils import *
from parsers import *

parser = argparse.ArgumentParser(description='Amalyzer of dataset')
parser.add_argument('interpro_set', help='File containing set of possible interpro ids')
parser.add_argument('ontology_file', help='File containing the go terms and their corresponding index') 
parser.add_argument('dataset', help='File containing the proteins to predict and their corresponding features')
parser.add_argument('-output_path', help='Path to the folder where the predction file will be put', default="InterPro_predictions")
args = parser.parse_args()

ip_set = args.interpro_set
go_file = args.ontology_file
data = args.dataset

ip_dict = parse_interpro_list(ip_set)

ontology = parse_ontology(go_file)

X , Y, ip_dict, go_dict, prot_dict = parse_dataset(data, ip_dict, ontology) 
X = X.toarray()
s = X.sum(axis=0)
#prots = prot_dict.keys()
#X = DataFrame(X, columns=list(ip_set))
#X["Protein id"] = prots

count = 0
for el in s:
    if el < 5 :
        count += 1
print("the number of interpro labels that appear less than 5 times is", count)

Y = Y.toarray()
sy = Y.sum(axis=0)
count = 0
for el in sy:
    if el < 5 :
        count += 1
print("the number of go terms that appear less than 5 times is", count)

rows, cols = X.shape
coupled_labels = []
for j in range(0,cols):
    l = list(range(0,cols))
    l.remove(j)
    coupled_labels.append(l)
for i in range(0,rows):
    for j in range(0,cols):
        for j2 in range(0,cols):
            if j != j2:
                if X[i,j] != X[i,j2]:
                    try:
                        coupled_labels[j].remove(j2)
                    except ValueError:
                        pass

print(coupled_labels)

                
        
        

