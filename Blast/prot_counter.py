from parsers import * 
from blast_score import *

def get_reference_pids(ref_file):
    """
    parses the ref file and returns a list of protein ids

    Parameters
    ----------
    ref_file: str 
        path to the file containing the reference 

    Returns
    -------
    uniprot : set
        set containing protein ids

    cafa : set
        set containing protein ids
    """
    uniprot = set()
    cafa = set()
    with open(ref_file, "r") as f:
        for line in f:
            div = line.split()
            uniprot.add(div[-1])
            cafa.add(div[0])
    return uniprot, cafa

def get_goa_pids(goa_file):
    """
    parses the ref file and returns a list of protein ids

    Parameters
    ----------
    goa_file: str 
        path to the file containing the goa db 

    Returns
    -------
    uniprot : set
        set containing protein ids
    """
    uniprot = set()
    with open(goa_file, "r") as f:
        for line in f:
            div = line.split()
            uniprot.add(div[0])
    return uniprot


def get_blast_pids(blast_file):
    """
    parses the ref file and returns a list of protein ids

    Parameters
    ----------
    goa_file: str 
        path to the file containing the goa db 

    Returns
    -------
    uniprot : set
        set containing protein ids
    """
    cafa = set()
    with open(blast_file, "r") as f:
        for line in f:
            div = line.split()
            cafa.add(div[0])
    return cafa


ref_uniprot, ref_cafa = get_reference_pids("/home/davide/Documenti/training/reference_30000.txt")
goa_uniprot = get_goa_pids("/home/davide/Documenti/training/goa_db_2018_08_exp.dat")
blast_cafa = get_blast_pids("/home/davide/Documenti/training/blast_pred_new.txt")

#common = ref_uniprot.intersection(goa_uniprot).intersection(blast_uniprot)
ref_goa_common = ref_uniprot.intersection(goa_uniprot)
ref_blast_common = ref_cafa.intersection(blast_cafa)
#ref_blast_common = ref_uniprot.intersection(blast_uniprot)
print("# of proteins in the reference:", len(ref_uniprot))
print("# of proteins in the goa file:", len(goa_uniprot))
#print("# of common ids between blast and ref:", len(ref_blast_common))
print("# of common ids between goa and ref:", len(ref_goa_common))
print("# of common ids between Blast and ref:", len(ref_blast_common))
#print("number of common protein:", len(common))