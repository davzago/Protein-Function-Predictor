import os
import argparse


parser = argparse.ArgumentParser(description='changes the predictions with uniprot ids to cafa ids')
parser.add_argument('predictions_folder', help='Path to the folder containing the predictions')
parser.add_argument('mapping', help="mapping between cafa ids and uniprot ids")
args = parser.parse_args()

preds = args.predictions_folder
mapping = args.mapping

id_map = dict()
with open(mapping, 'r')as f:
        for line in f:
            cafa_id, uniprot_id , _ = line.split()
            id_map.setdefault(uniprot_id, cafa_id)

for pred in os.listdir(preds):
    pred_list = []
    with open(preds + '/' + pred, 'r') as f:
        for line in f:
            uniprot_id, go_term , score = line.split()
            pred_list.append([uniprot_id, go_term, score])
    with open(preds + '/' + pred, 'w') as f:
        for uniprot_id, go_term, score in pred_list:
            if uniprot_id in id_map:
                cafa_id = id_map[uniprot_id] 
                f.write(cafa_id + '\t' + go_term + '\t' + score + '\n')
            