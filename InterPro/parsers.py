import time

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

def parse_features(dataset_file, interpro_set, ontology):
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
    x_train : csr_matrix
        a compressed sparse row matrix representing the features for each protein

    y_train : csr_matrix
        a compressed sparse row matrix representing the ground truth for each go term

    ip_dict : dict
         {interpro_id : index}

    prot_dict : dict
        {uniprot_id : (cafa_id, index)}
    """
    ip_dict = set_ip_indices(interpro_set)
    go_dict = set_ont_indices(ontology)
    prot_indexes = dict()
    idx = 0
    x_row = []
    x_col = []
    x_data = []
    y_row = []
    y_col = []
    y_data = []

    X = []

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
            # fill the ground truth matrix
            for go_id in go_id_list.split('-'):
                y_row.append(idx)
                y_col.append(go_dict[go_id])
                y_data.append(1)
            
            idx += 1
    
    x_row = np.array(x_row)
    x_col = np.array(x_col)
    x_data = np.array(x_data, dtype='bool')
    y_row = np.array(y_row)
    y_col = np.array(y_col)
    y_data = np.array(y_data, dtype='bool')

    n_prot = idx
    n_ip_id = len(interpro_set)
    n_go_terms = len(ontology) 

    x_train = csr_matrix((x_data, (x_row, x_col)), (n_prot, n_ip_id))
    y_train = csc_matrix((y_data, (y_row, y_col)), (n_prot, n_go_terms))

    return x_train, y_train, ip_dict, go_dict, prot_indexes