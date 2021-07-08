import numpy as np
from scipy.sparse import csr_matrix

def create_dataset(interPro_file, reference_dict, ontology, interpro_ids_set):
    """
    Takes interpro assocs and a reference file and returns a dataset  
    NB since n = #go terms classifiers will be needed the dataset will contain n label cols

    Parameters
    ----------
    interPro_file : str
        path to the interpro file
    
    reference_dict : dict
        {uniprot_id: (cafa_id, set(go_terms))}

    ontology : dict
        {go_id : list(parents)}

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

    go_dict : dict
        {go_id : index}
    """
    x_row = []
    x_col = []
    x_data = []
    y_row = []
    y_col = []
    y_data = []
    prot_idx = 0
    prot_indexes = dict()

    ip_dict = set_ip_indices(interpro_ids_set)
    go_dict = set_ont_indices(ontology)
    # cycle on interpro file and check of is are contained in the reference
    # if so add the protein to the dataset by adding row and col indices to the arrays declared above
    with open(interPro_file) as f:
        for line in f:
            uniprot_id, interpro_list = line.split()
            if uniprot_id in reference_dict:
                prot_indexes.setdefault(uniprot_id, prot_idx)
                # fill up matrix containing the features x_train
                for ip_id in interpro_list.split('-'):
                    x_row.append(prot_idx)
                    x_col.append(ip_dict[ip_id])
                    x_data.append(1)
                # fill up the ground truth matrix y_train
                for go_id in reference_dict[uniprot_id][1]:
                    if go_id in ontology:
                        y_row.append(prot_idx)
                        y_col.append(go_dict[go_id])
                        y_data.append(1)
                prot_idx += 1
    
    # create the actual compressed matrices
    x_row = np.array(x_row)
    x_col = np.array(x_col)
    x_data = np.array(x_data, dtype='bool')
    y_row = np.array(y_row)
    y_col = np.array(y_col)
    y_data = np.array(y_data, dtype='bool')

    n_prot = prot_idx
    n_ip_id = len(interpro_ids_set)
    n_go_terms = len(ontology) 

    x_train = csr_matrix((x_data, (x_row, x_col)), (n_prot, n_ip_id))
    y_train = csr_matrix((y_data, (y_row, y_col)), (n_prot, n_go_terms))

    return x_train, y_train, ip_dict, go_dict

    
    
def parse_reference(ref_file):
    """
    Takes a reference file, parses it and returns the corresponding dictionary

    Parameter
    ---------
    ref_file : str
        Path to the reference file 

    Returns
    -------
    ref_dict : dict
        {uniprot_id: (cafa_id, set(go_terms))}
    """
    ref_dict = dict()
    with open(ref_file) as f:
        for line in f:
            cafa_id, go_id, namespace, uniprot_id = line.split()
            ref_dict.setdefault(uniprot_id, (cafa_id, set()))
            ref_dict[uniprot_id][1].add(go_id)
    return ref_dict

def set_ip_indices(interpro_set):
    """
    Takes a set containing all the possible intepro ids and returns a dictionari that assgins to each id a index

    Parameters
    ----------
    interpro_set : set
        set containing all the interpro_terms

    Returns
    -------
    index_dict : dict
        {interpro_id : index}
    """
    index_dict = dict()
    index = 0
    for ip_id in interpro_set:
        index_dict.setdefault(ip_id, index)
        index += 1
    return index_dict

def set_ont_indices(ontology):
    """
    Takes a set containing all the possible intepro ids and returns a dictionari that assgins to each id a index

    Parameters
    ----------
    ontology : dict
        {go_id : list(parents)}

    Returns
    -------
    index_dict : dict
        {go_id : index}
    """
    index_dict = dict()
    index = 0
    for go_id in ontology.keys():
        index_dict.setdefault(go_id, index)
        index += 1
    return index_dict

def reverse_dict(dictionary):
    """
    Takes a dictionary where the values are indices and returns a list where each index position contains the key

    Parameters
    ----------
    dictionary : dict
        {id : index}

    Returns
    -------
    rev : np.array
        array containing for each index the corresponding id
    """
    rev = np.full(len(dictionary), '')
    for k, v in dictionary.items():
        rev[v] = k

    return rev

def save_dataset(x_train, y_train, ip_dict, go_dict, output_path):
    """
    takes the dataset and the corresponding dicts and saves it on a text file

    Parameters
    ----------
    x_train : csr_matrix
        matrix with rows representing proteins and columns representing interpro ids

    y_train : csr_matrix
        matrix with rows representing proteins and columns representing go_terms

    ip_dict : dict
        {ip_id : index}

    go_dict : dict
        {go_id : index}

    output_path : str
        Path where to put the output fileg
    """



