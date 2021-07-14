import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix

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

    prot_dict : dict
        {uniprot_id : index}
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
    y_train = csc_matrix((y_data, (y_row, y_col)), (n_prot, n_go_terms))

    return x_train, y_train, ip_dict, go_dict, prot_indexes

    
    
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
        {uniprot_id: [cafa_id, set(go_terms)]}
    """
    ref_dict = dict()
    with open(ref_file) as f:
        for line in f:
            cafa_id, go_id, namespace, uniprot_id = line.split()
            ref_dict.setdefault(uniprot_id, [cafa_id, set()])
            ref_dict[uniprot_id][1].add(go_id)
    return ref_dict

def propagate_reference(ref_dict, ont):
    """
    Takes a reference dictionary and propagates the go terms associated with each protein

    Parameters
    ----------
    ref_dict : dict
        {uniprot_id: (cafa_id, set(go_terms))}

    ontology : dict
        {go_id : list(parents)}

    Returns
    -------
    ref_dict : dict
        {uniprot_id: (cafa_id, set(go_terms))}
    """
    for k, v in ref_dict.items():
        v[1] = propagate_terms(v[1], ont)

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
    rev = np.empty(len(dictionary), dtype=object)
    for k, v in dictionary.items():
        rev[v] = k 

    return rev

def save_dataset(x_train, y_train, prot_dict, ip_dict, go_dict, output_path, ref_dict):
    """
    takes the dataset and the corresponding dicts and saves it on a text file

    Parameters
    ----------
    x_train : csr_matrix
        matrix with rows representing proteins and columns representing interpro ids

    y_train : csr_matrix
        matrix with rows representing proteins and columns representing go_terms

    prot_dict : dict
        {prot_id : index}

    ip_dict : dict
        {ip_id : index}

    go_dict : dict
        {go_id : index}

    output_path : str
        Path where to put the output file
    """

    rev_ip = reverse_dict(ip_dict)
    rev_go = reverse_dict(go_dict)
    rev_prot = reverse_dict(prot_dict)


    n, m  = x_train.shape
    p = y_train.shape[1]
    f = open(output_path+"/dataset.txt", "w")
    for i in range(0,n):
        f.write(rev_prot[i] + "\t" + ref_dict[rev_prot[i]][0] + "\t")
        indexes = x_train[i,:].nonzero()
        pos = -1 
        for idx in indexes[1]:
            pos += 1
            ip_id = rev_ip[idx]
            f.write(ip_id)
            if pos != len(indexes[1]) - 1:
                f.write('-')
        f.write('\t')
        indexes = y_train[i,:].nonzero()
        pos = -1
        for idx in indexes[1]:
            pos += 1
            go_id = rev_go[idx]
            f.write(go_id)
            if pos != len(indexes[1]) - 1:
                f.write('-')
        if i != n-1:
            f.write("\n")
    f.close()
    

def parse_dataset(dataset_file, interpro_set, ontology):
    """
    Parses a dataset file created with the function save_dataset() and returns the same structures as create_dataset()

    Parameters
    ----------
    dataset_file : str
        path to the dataset file

    interpro_ids_set : set
        set containing all the interpro ids

    ontology : dict
        {go_id : list(parents)}
    
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


def save_interpro_set(interpro_set, output_path):
    """
    Saves on a text file a set containing all the possible interpro ids

    Parameters
    ----------
    interpro_set : set
        set of interpro ids

    output_path : str
        path where to put the output text file 
    """
    f = open(output_path + "/interpro_set.txt", "w")
    for ip_id in interpro_set:
        f.write(ip_id + '\n')
    f.close()

def parse_interpro_set(ip_set_file):
    """
    Parses a text file containing all the possible interpro ids and returns the corresponding set

    Parameters
    ----------
    ip_set_file : ste
        path to the interpro ids file

    Returns
    -------
    ip_set : set
        set containing all the possible interpro ids
    """
    ip_set = set()
    with open(ip_set_file) as f:
        for line in f:
            ip_set.add(line.strip('\n'))
    return ip_set

def propagate_terms(go_term_set, ont):
    """
    Takes a go term set and adds all the ancestors of each term

    Parameters
    ----------
    go_term_set : set
        set of go terms
    
    ont : dict
        {go_id: list(parents)}

    Returns
    -------
    go_term_set
        set containing go terms and ancestors of the initial input
    """
    queue = list(go_term_set)
    while queue:
        term = queue.pop(0)
        if term in ont:
            parents = ont[term]
            for p in parents:
                go_term_set.add(p)
                queue.append(p)
    return go_term_set

def remove_unused_go_terms(y_train, y_test, go_dict):
    """
    removes from the ground truth the columns that ar all zero

    Parameters
    ----------
    y_train : matrix
        matrix where rows correspond to proteins and coumns correspond to go terms

    y_test : matrix
        matrix where rows correspond to proteins and coumns correspond to go terms

    go_dict : dictionary
        {go_term: index}

    Returns
    -------
    y : matrix
        input matrix without 0 cols
        
    go_dict : dictionary
        {go_term: index} with rearranged index
    """

    go_arr = reverse_dict(go_dict)
    cols = list(set(np.nonzero(y_train)[1]))
    y_new = csc_matrix(y_train[:,cols])
    test = csc_matrix(y_test[:,cols])
    go_arr_new = go_arr[cols]

    go_dict_new = dict()
    for i in range(0,len(go_arr_new)):
        go_dict_new[go_arr_new[i]] = i

    return y_new, test, go_dict_new


def save_dict(dict, output_path, file_name):
    """
    Takes a dictionary and saves it on a text file 

    Parameters
    ----------
    dict : dictionary
        {id : index}

    output_path : str
        path to the output folder

    file_name : str
        name of the text file

    """
    with open(output_path + "/" + file_name + ".txt", "w") as f:
        for k, v in dict.items():
            f.write(k + '\t' + str(v) + '\n')
            

def save_prot_dict(dict, output_path, file_name):
    """
    Takes a dictionary of proteins and saves it on a text file 

    Parameters
    ----------
    dict : dictionary
        {id : index}

    output_path : str
        path to the output folder

    file_name : str
        name of the text file

    """
    with open(output_path + "/" + file_name + ".txt", "w") as f:
        for k, v in dict.items():
            f.write(k + '\t' + v[0] + '\t' + str(v[1]) + '\n')





    