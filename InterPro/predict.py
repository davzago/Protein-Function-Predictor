from InterPro.parsers import parse_ontology
from InterPro.dataset_utils import parse_interpro_set
import argparse
from joblib import load



parser = argparse.ArgumentParser(description='Predictor based on InterPro')
parser.add_argument('interpro_set', help='File containing set of possible interpro ids')
parser.add_argument('ontology_file', help='File containing the ontology')
parser.add_argument('features_file', help='File containing the the interPro feature of the of each protein')
parser.add_argument('-output_path', help='Path to the folder where the predction file will be put', default="output")
args = parser.parse_args()

interpro_file = args.interpro_file
features_file = args.features_file
ontology_file = args.ontology_file
output_path = args.output_path

ip_set = parse_interpro_set(interpro_file)
ontology = parse_ontology(ontology_file)
#X = parse_features()

