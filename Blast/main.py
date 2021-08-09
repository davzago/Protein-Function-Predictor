from parsers import *
from blast_score import *
import argparse
import os


# If a uniprot id is in both reference and goa then the results will be rigged because basically goa will say what go term to assign

parser = argparse.ArgumentParser(description='Blast predictor')
parser.add_argument('ontology_file', help='File containing the ontology')
parser.add_argument('goa_file', help='File containing the goa database')
parser.add_argument('blast_file', help="File containing the blast results")
parser.add_argument('reference_file', help="File that contains all the proteins that have to be predicted")
parser.add_argument('-output_path', help='Path to the folder where the predction file will be put', default="output")
args = parser.parse_args()


goa_file = args.goa_file
output_path = args.output_path
blast_file = args.blast_file
ref_file = args.reference_file
ontology_file = args.ontology_file

if not os.path.isdir(output_path):
    os.mkdir(output_path)

ref_set = get_query_ids(ref_file)
ont_dict, _ = parse_ontology(ontology_file)
present_dict = parse_blast(blast_file, goa_file, ont_dict, ref_set, output_path)
present_dict = reverse_dict(present_dict)
present_dict = parse_goa(goa_file, present_dict, ont_dict)
present_dict, gt_dict, go_term_set = return_to_query(present_dict)
present_dict = blast_score_no_out(present_dict, gt_dict, go_term_set)

save_predictions(present_dict, "predictions", output_path)
