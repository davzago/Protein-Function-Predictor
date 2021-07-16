import time
from dataset_utils import *

def interpro_ids(IPro_file):
    """
    Takes an already parsed InterProScan result file and returns the total number of unique accessions

    Parameters
    ----------
    IPro_file : str
        file containing the parsed results of an InterProScan

    Returns
    -------
    set_length : set
        set containing all the possible interpro ids
    """
    ip_set = set()
    with open(IPro_file) as f:
        for line in f:
            uniprot_id, accessions = line.split()
            acc = accessions.split('-')
            for e in acc:
                ip_set.add(e)
    return ip_set

"""start = time.time()
unique_ipro = count_unique_accessions("/home/davide/Documenti/training/interpro_parsed.txt")
print("the number of unique ids is:", unique_ipro)
print("time elapsed", time.time()-start)"""

def parse_ontology(obo_file_path):
    """
    Parses an ontology file and returns the corresponding dictionary

    Parameters
    ----------
    obo_file_path : str
        path to the ontology file:

    Returns
    -------
    go_dict : dict
        {go_id : list(parents)}
    """
    go_terms = dict()
    alts = dict()
    temp_id = ''
    temp_alt_id = []
    rel_list = []
    with open(obo_file_path) as f:
        rel_lost =  []
        for line in f:
            line = line.strip().split(": ")
            if line and len(line) > 1:
                k, v = line[:2]
                if k == "id" and (v.startswith("GO:") or v.startswith("HP:") or v.startswith("DO:")):
                    # when a new id is found first we have to input the entry in the go term dict
                    # also we will have to input a new entry when we reach EOF 
                    if temp_id != '':
                        go_terms[temp_id] = rel_list
                        alts[temp_id] = temp_alt_id
                        rel_list = []
                        temp_alt_id = []
                    temp_id = v
                elif k == "alt_id" and (v.startswith("GO:") or v.startswith("HP:") or v.startswith("DO:")):
                    temp_alt_id.append(v)
                elif k == "is_a" and (v.startswith("GO:") or v.startswith("HP:") or v.startswith("DO:")):
                    # add (temp_id,v) tuple to the relation list
                    s = v.split('!')[0].strip()
                    rel_list.append(s)

        go_terms[temp_id] = rel_list
        alts[temp_id] = temp_alt_id

    return go_terms

def parse_features(dataset_file, interpro_set):
    """
    Parses a dataset file created with the function save_dataset() and returns the same structures as create_dataset()

    Parameters
    ----------
    dataset_file : str
        path to the dataset file

    interpro_ids_set : set
        set containing all the interpro ids

    Returns
    -------
    X : csr_matrix
        a compressed sparse row matrix representing the features for each protein

    ip_dict : dict
         {interpro_id : index}

    prot_dict : dict
        {uniprot_id : (cafa_id, index)}
    """
    
    ip_dict = set_ip_indices(interpro_set)
    prot_indexes = dict()
    idx = 0
    x_row = []
    x_col = []
    x_data = []

    with open(dataset_file) as f:
        for line in f:
            #rowx = np.zeros(shape)
            uniprot_id, cafa_id, interpro_list, go_id_list = line.split()
            prot_indexes.setdefault(uniprot_id, (cafa_id, idx))
            # fill the training set matrix
            for ip_id in interpro_list.split('-'):
                x_row.append(idx)
                x_col.append(ip_dict[ip_id])
                x_data.append(1)
            idx += 1
    
    x_row = np.array(x_row)
    x_col = np.array(x_col)
    x_data = np.array(x_data, dtype='bool')

    n_prot = idx
    n_ip_id = len(interpro_set)

    x_train = csr_matrix((x_data, (x_row, x_col)), (n_prot, n_ip_id))

    return x_train, ip_dict, prot_indexes

def parse_dict(file_path):
    idx_dict = dict()
    with open(file_path) as f:
        for line in f:
            k, v = line.split()
            idx_dict[k] = v
    return idx_dict

def parse_prot_dict(file_path):
    prot_dict = dict()
    with open(file_path) as f:
        for line in f:
            uniprot_id, cafa_id, index = line.split()
            prot_dict[uniprot_id] = (cafa_id, index)
    return prot_dict

def parse_dict_to_array(file_path):
    idx_arr = []
    with open(file_path) as f:
        for line in f:
            k, v = line.split()
            idx_arr.append(k)
    idx_arr = np.array(idx_arr, dtype=object)
    return idx_arr
