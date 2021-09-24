from parsers import *
from blast_score import *
import argparse
import os
import pickle


# If a uniprot id is in both reference and goa then the results will be rigged because basically goa will say what go term to assign

"""parser = argparse.ArgumentParser(description='Blast predictor')
parser.add_argument('goa_file', help='File containing the goa database')
parser.add_argument('reference_file', help="File that contains all the matches between cafa id and uniprot id")
parser.add_argument('-output_path', help='Path to the folder where the predction file will be put', default="output")
args = parser.parse_args()


goa_file = args.goa_file
output_path = args.output_path
ref_file = args.reference_file"""

def connect_ids(references):
    match_dict = dict()
    rev_match = dict()
    for r in references:
        with open(r, 'r') as f:
            for line in f:
                cafa_id, _, _, uniprot_id = line.split() 
                match_dict.setdefault(uniprot_id, cafa_id)
                rev_match.setdefault(cafa_id, uniprot_id)
    return match_dict, rev_match

def match_uniprot_ids(goa, match_dict):
    ids = set()
    count = 0
    with open(goa, 'r') as f:
        for line in f:
            uniprot_id = line.split()[0]
            if uniprot_id in match_dict:
                ids.add(match_dict[uniprot_id])
            count += 1
    print("# of proteins:", count)
    return ids


refs = ["/home/davide/Documenti/training/reference_30000.txt", "/home/davide/Documenti/training/reference_new.txt"]
match_dict, rev_match = connect_ids(refs)
ids = match_uniprot_ids("/home/davide/Documenti/training/goa_db_2018_08_exp.dat", match_dict)
print(len(ids))

dict_file =open("ids_dict.pkl", "wb")

pickle.dump(match_dict, dict_file)

