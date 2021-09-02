from parsers import * 
from blast_score import *
from utils import *


ref_uniprot, ref_cafa = get_reference_pids("/home/davide/Documenti/training/reference_new.txt")
goa_uniprot = get_goa_pids("/home/davide/Documenti/training/goa_db_2018_08_exp.dat")
blast_cafa = get_blast_pids("/home/davide/Documenti/training/blast_pred.txt")

#common = ref_uniprot.intersection(goa_uniprot).intersection(blast_uniprot)
ref_goa_common = ref_uniprot.intersection(goa_uniprot)
ref_blast_common = ref_cafa.intersection(blast_cafa)
#ref_blast_common = ref_uniprot.intersection(blast_uniprot)
print("# of proteins in the reference:", len(ref_uniprot))
print("# of proteins in the goa file:", len(goa_uniprot))
#print("# of common ids between blast and ref:", len(ref_blast_common))
print("# of common ids between goa and ref:", len(ref_goa_common))
print(list(ref_goa_common)[0])
print("# of common ids between Blast and ref:", len(ref_blast_common))
#print("number of common protein:", len(common))