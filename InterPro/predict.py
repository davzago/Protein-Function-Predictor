from parsers import parse_dict, parse_dict_to_array
from dataset_utils import parse_interpro_set, reverse_dict, save_prediction, set_ip_indices
import argparse
import numpy as np

from joblib import load



parser = argparse.ArgumentParser(description='Predictor based on InterPro')
parser.add_argument('interpro_set', help='File containing set of possible interpro ids')
parser.add_argument('ontology_dict_file', help='File containing the go terms and their corresponding index') 
parser.add_argument('model', help='Path to the model used to output the predictions')
parser.add_argument('protein_features', help='String containing the uniprot id, cafa id and associated interpro ids')
parser.add_argument('-output_path', help='Path to the folder where the predction file will be put', default="InterPro_predictions")
args = parser.parse_args()

interpro_set = args.interpro_set
features = args.protein_features
terms_dict = args.ontology_dict_file
model = args.model
output_path = args.output_path

clf = load(model)
ip_set = parse_interpro_set(interpro_set)
ip_dict = set_ip_indices(ip_set)
#go_dict = parse_dict(terms_dict)
#go_array = reverse_dict(go_dict)
go_array = parse_dict_to_array(terms_dict)

uniprot_id, cafa_id, ip_list = features.split()
x = np.zeros(len(ip_dict), dtype='bool')

for ip_id in ip_list.split('-'):
    x[ip_dict[ip_id]] = 1


pred = clf.predict_proba(go_array)

go_scores = dict()
for idx, score in enumerate(pred):
    if score > 0.01:
        go_term = go_array[idx]
        go_scores[go_term] = score

save_prediction(cafa_id, go_scores, output_path)







