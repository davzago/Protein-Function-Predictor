from numpy.lib.npyio import save
from dataset_utils import *
from parsers import *


interpro_file = "/home/davide/Documenti/training/interpro_parsed.txt"
ref_file = "/home/davide/Documenti/training/reference_new.txt"
ontology_file = "/home/davide/Documenti/training/gene_ontology_edit.obo.2018-08-01"

ip_set = interpro_ids(interpro_file)
save_interpro_set(ip_set, "InterPro")
ontology = parse_ontology(ontology_file)
ref_dict = parse_reference(ref_file)

x_train, y_train, ip_dict, go_dict, prot_dict = create_dataset(interpro_file, ref_dict, ontology, ip_set)

save_dataset(x_train, y_train, prot_dict, ip_dict, go_dict, "InterPro", ref_dict)