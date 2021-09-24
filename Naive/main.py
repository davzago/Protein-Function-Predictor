import argparse
import os
from naive_functions import *



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

if not os.path.isdir(output_path):
    os.mkdir(output_path)

ontology = parse_ontology(obo_file)
goa_dict = parse_goa(goa_file)
freq_dict = naive_calculator(goa_dict, ontology)
#prot_set = cafaid_form_ref(protein_file)
prot_set = get_protein_set(protein_file)
naive_to_prediction(freq_dict, prot_set, output_path)





