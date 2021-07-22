from numpy.lib.npyio import save
from dataset_utils import *
from parsers import *
import argparse
import os

parser = argparse.ArgumentParser(description='dataset builder')
parser.add_argument('ontology_file', help='File containing the ontology')
parser.add_argument('interpro_file', help='File containing the goa database')
parser.add_argument('ref_file', help="File containing the cafa id of the proteins that you want to predict")
parser.add_argument('-output_gt', help="if true returns a dataset file containing features and gt, if false returns a file containing a feature matrix",
                    default=True)
parser.add_argument('-output_path', help='Path to the folder where the predction file will be put', default="InterPro")
args = parser.parse_args()

interpro_file = args.interpro_file
ref_file = args.ref_file
ontology_file = args.ontology_file
output_path = args.output_path
gt = args.output_gt

if not os.path.isdir(output_path):
    os.mkdir(output_path)

ip_set = interpro_ids(interpro_file)
save_interpro_set(ip_set, output_path)
ontology = parse_ontology(ontology_file)
ref_dict = parse_reference(ref_file)
ref_dict = propagate_reference(ref_dict, ontology)

x_train, y_train, ip_dict, go_dict, prot_dict = create_dataset(interpro_file, ref_dict, ontology, ip_set)

if gt == True:
    save_dataset(x_train, y_train, prot_dict, ip_dict, go_dict, output_path, ref_dict)
else:
    save_reference(x_train, prot_dict, ip_dict, output_path, ref_dict)