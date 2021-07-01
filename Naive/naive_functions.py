

def parse_ontology(obo_file_path):
    """
    Parses a file containing an ontology and returns the corresponding dictionary

    Parameters
    ----------
    obo_file_path : str
        the location of the ontology file

    Returns
    -------
    dictionary
        {go_id: list(parents)}
    """
    go_terms = dict()
    alts = dict()
    temp_id = ''
    temp_alt_id = []
    rel_list = []
    with open(obo_file_path) as f:
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

    return go_terms#, alts



def parse_goa(goa_file):
    """
    Parses a goa file and returns its data as a dictionary

    Parameters
    ----------
    goa_file : str
        the file location of the goa_db

    Returns
    -------
    dictionary
        {uniprot_id: set(go_id)}
    """
    goa_dict = dict()
    with open(goa_file) as f:
        for line in f:
            uniprot_id, _, terms = line.split()
            for go_list in terms.split(';'):
                l = go_list.split(',')
                for i in range(1,len(l)):
                    goa_dict.setdefault(uniprot_id, set())
                    goa_dict[uniprot_id].add(l[i])
    return goa_dict

def propagate_terms(go_terms_set, ontology):
    """
    takes a set of go terms and propagates adds the ancestors of each 

    Parameters
    ----------
    go_terms : set
        set of go terms

    ontology : dict
        dictionary containing the ontology

    Returns
    -------
    set
        set of propagated go terms
    """
    queue = list(go_terms_set)
    while queue:
        term = queue.pop(0)
        if term in ontology:
            parents = ontology[term]
            for p in parents:
                go_terms_set.add(p)
                queue.append(p)
    return go_terms_set


def naive_calculator(goa_dict, ontology):
    """
    Computes the frequency of go terms in a dataset

    Parameters
    ----------
    goa_dict : dictionary
        {uniprot_id: set(go_id)}

    ontology : dictionary
        {go_id: list(parents)}

    Returns
    -------
    dictionary
        {go_id: frequency}
    """
    frequency_dict = dict()
    for go_set in goa_dict.values():
        go_set = propagate_terms(go_set, ontology)
        for go_id in go_set:
            frequency_dict.setdefault(go_id, 0)
            frequency_dict[go_id] += 1
    
    N_d = len(goa_dict)
    for go_id in frequency_dict.keys():
        frequency_dict[go_id] /= N_d
    
    return frequency_dict

def naive_to_prediction(frequency_dict, protein_list, out_path, max_n_prot=500):
    """
    Transforms the friquency of go_terms in predictions ready for evaluation

    Parameters
    ----------
    frequency_dict : dictionary
        {go_id: frequency}

    protein_list : list
        list of cafa ids
    
    output_path : str
        path to the output folder

    max_n_prot : int
        maximum number of proteins per prediction file
    """
    i = 1
    count = 0
    f = open(out_path+ '/' + "predictions" + str(i) + ".txt", "w")
    for p_id in protein_list:
        count += 1
        if count > max_n_prot:
            i = i+1
            f.close()
            f = open(out_path+ '/' + "predictions" + str(i) + ".txt", "w")
            count = 0
        for go_id, freq in frequency_dict.items():
            if freq > 0.01:
                f.write(p_id + "\t" + go_id + "\t" + str(round(freq, 2)) + '\n')
    f.close()

def cafaid_form_ref(file):
    """
    returns all the protein cafa id in a reference file

    Parameters
    ----------
    file : str
        path to a file containing the reference

    Returns
    -------
    set
        a set containing all the protein cafa ids in the reference file
    """
    cafa_ids = set()
    with open(file) as f:
        for line in f:
            cafa_id = line.split()[0]
            cafa_ids.add(cafa_id)
    return cafa_ids

"""def benchmark_from_ref(file, gt_file):
    cafa_ids = set()
    gt_ids = set()
    with open(file) as f:
        for line in f:
            cafa_id = line.split()[0]
            cafa_ids.add(cafa_id)
    
    with open(gt_file) as gt:
        for line in gt:
            c_id = line.split()[0]
            gt_ids.add(c_id)


    wf = open("benchmark.txt", "w")
    for cafa_id in cafa_ids:
        if cafa_id in gt_ids:
            wf.write(cafa_id + '\n')
    wf.close()"""

#cafa_ids = cafaid_form_ref("/home/davide/Documenti/training/reference_new.txt")
#benchmark_from_ref("/home/davide/Documenti/training/reference_new.txt", "/home/davide/Documenti/Protein Function Predictor/Blast/gt_folder/CCO_new_gt.txt")



        

        
