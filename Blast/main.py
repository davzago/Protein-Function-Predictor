from parsers import *
from blast_score import *
import argparse
import os

parser = argparse.ArgumentParser(description='Blast predictor')
parser.add_argument('ontology_file', help='File containing the ontology')
parser.add_argument('goa_file', help='File containing the goa database')
parser.add_argument('blast_file', help="File containing the blast results")
parser.add_argument('-output_path', help='Path to the folder where the predction file will be put', default="output")
args = parser.parse_args()


goa_file = args.goa_file
output_path = args.output_path
blast_file = args.blast_file
ontology_file = args.ontology_file

if not os.path.isdir(output_path):
    os.mkdir(output_path)


ont_dict, _ = parse_ontology(ontology_file)
blast_dict = parse_blast(blast_file, goa_file, ont_dict, output_path)