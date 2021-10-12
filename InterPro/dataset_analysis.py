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

print(max(s))
for key, idx in ip_dict.items():
    if idx == 0:
        print(key)
