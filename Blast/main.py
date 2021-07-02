from parsers import *
from blast_score import *
import argparse

parser = argparse.ArgumentParser(description='Naive predictor')
parser.add_argument('ontology_file', help='File containing the ontology')
parser.add_argument('goa_file', help='File containing the goa database')
parser.add_argument('protein_file', help="File containing the cafa id of the proteins that you want to predict")
parser.add_argument('-output_path', help='Path to the folder where the predction file will be put', default="output")
args = parser.parse_args()


goa_file = args.goa_file
output_path = args.output_path
protein_file = args.protein_file
obo_file = args.ontology_file


ont_dict, _ = parse_ontology(obo_file)
blast_dict = parse_blast("fake_blast.txt")
"""print("blast parsed")
rev = reverse_dict(blast_dict)
blast_dict.clear()
print("dictionary reversed")
goa_dict = parse_goa("/home/davide/Documenti/training/goa_db_2018_08_exp.dat", rev, ont_dict)
rev.clear()
print("goa parsed")
query_dict, gt_dict, go_term_set = return_to_query(goa_dict)
goa_dict.clear()
print("back to query")
score_dict = blast_score(query_dict, gt_dict, go_term_set)
print("score calculation")
prediction_to_text(score_dict, 100)"""
#parse_reference("/home/davide/Documenti/training/reference_new.txt")