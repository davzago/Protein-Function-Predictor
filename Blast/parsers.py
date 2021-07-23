from blast_score import *
import gzip

# could include part_of of the parents
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

    return go_terms, alts

# takes a reference file and splits it in a file for each namespace
# (the uniport ids are also removed because they aren't useful to the evaluation)
def parse_reference(reference_file):
    """
    parses a reference file and returns a file for each different namespace

    Parameters
    ----------
    reference_file : str
        path to the reference file
    """
    namespaces = dict()
    with open(reference_file) as rf:
        for line in rf:
            cafa_id, go_id, namespace, _ = line.split()
            if namespace in namespaces:
                namespaces[namespace].append([cafa_id, go_id])
            else:
                namespaces.setdefault(namespace, [[cafa_id, go_id]])
    rf.close()
    for ns, ref_list in namespaces.items():
        f = open(ns+"_gt.txt", "w")
        for cafa_id, go_id in ref_list:
            f.write(cafa_id + "\t" + go_id + "\n")
        f.close()


def parse_blast(blast_file, goa_file, ont_dict, output_path):
    """
    Parses the blast file in chunks and return and writes the prediction files using blast_score

    Parameters
    ----------
    blast_file : str
        path to the blast file

    goa_file : str
        path to the goa file

    ont_dict : dict
        {go_id: list(parents)}

    output_path : str
        path to the output folder 
    """
    blast_dict = dict()
    with open(blast_file) as f:
        prev_cafa_id = ''
        count = 0
        i = 1
        for line in f:
            cafa_id, uniprot_id, _, _, _, _, _, _, _, _, e_value, bit_score = line.split()
            ### When having memory problems uncomment ths code and run only this function (will write the prediction file)
            if prev_cafa_id == '':
                prev_cafa_id = cafa_id
            if prev_cafa_id != cafa_id:
                count += 1
                if count > 1000:
                    rev = reverse_dict(blast_dict)
                    goa_dict = parse_goa(goa_file, rev, ont_dict)
                    query_dict, gt_dict, go_term_set = return_to_query(goa_dict)
                    blast_score(query_dict, gt_dict, go_term_set, output_path, i)
                    i += 1
                    count = 0
                prev_cafa_id = cafa_id
            if float(e_value) < 0.001:
                if cafa_id in blast_dict:
                    blast_dict[cafa_id].append([uniprot_id, float(e_value), float(bit_score)])
                else:
                    blast_dict.setdefault(cafa_id, [[uniprot_id, float(e_value), float(bit_score)]])
        #This aswell
        rev = reverse_dict(blast_dict)
        goa_dict = parse_goa(goa_file, rev, ont_dict)
        query_dict, gt_dict, go_term_set = return_to_query(goa_dict)
        blast_score(query_dict, gt_dict, go_term_set, output_path, i)


def reverse_dict(blast_dict):
    """
    takes a dictionary with cafa_id as key and reverses it so the key is uniprot_id

    Parameters
    ----------
    blast_dict : dict
        {cafa_id: [uniprot_id, e_value, bit_score]}

    Returns
    -------
    rev_dict : dict
        {uniprot_id: (go_terms_set,[cafa_id, e_value, bit_score])}
    """
    rev_dict = dict()
    for cafa_id in list(blast_dict):
        for uniprot_id, e_value, bit_score in blast_dict[cafa_id]:
            rev_dict.setdefault(uniprot_id, [set(),[]])
            rev_dict[uniprot_id][1].append([cafa_id, e_value, bit_score])
        del blast_dict[cafa_id]
    return rev_dict

#Parses a GOA file and adds to the uniprot_dict(obtained in reverse_dict) a set of associated go-terms
def parse_goa(goa_file, rev_dict, ont_dict):
    """
    Parses a GOA file and adds to the uniprot_dict(obtained in reverse_dict) a set of associated go-terms

    Parameters
    ----------
    goa_file : str
        Path to the file containing the goa db

    rev_dict : dict
        {uniprot_id: (go_terms_set,[cafa_id, e_value, bit_score])}

    ont_dict : dict
        {go_id : list(parents)}

    Returns
    -------
    rev_dict : dict
        same dictionary taken in input with the go terms added and propagated in the set

    """
    with open(goa_file) as f:
        for line in f:
            uniprot_id, _, terms = line.split()
            for go_list in terms.split(';'):
                l = go_list.split(',')
                for i in range(1,len(l)):
                    if uniprot_id in rev_dict:
                        rev_dict[uniprot_id][0].add(l[i])
                        #print(uniprot_id)
                    #else:
                        #print(uniprot_id, "is not contained in the dictionary")
    k_list = []
    for k, v in rev_dict.items():
        if len(v[0]) == 0:
            k_list.append(k)
        else:
            v[0] = propagate_terms(v[0], ont_dict) # PROPAGATION
            
    for k in k_list:
        rev_dict.pop(k)
    return rev_dict

# takes a dictionary uni_dict = {uniprot_id: (go_terms_set,[cafa_id, e_value, bit_score])} and
# returns a dictionary query_dict = {query_id: {uniprot_id: bit_score}} and a gt_dict = {uniprot_id: go_term_set}
def return_to_query(sbj_dict):
    """
    takes a dictionary with uniprot id as key and switches it back to query id(cafa id) as key

    Parameters
    ----------
    sbj_dict : dict
        {uniprot_id: (go_terms_set,[cafa_id, e_value, bit_score])}

    Returns
    -------
    query_dict : dict 
        {query_id: {uniprot_id: bit_score}}

    Gt_dict : dict
        {uniprot_id: go_term_set}

    go_term_set : set
        set containing all the go terms found in the goa database   
    """
    query_dict = dict()
    gt_dict = dict()
    go_term_set = set()
    for k, v in sbj_dict.items():
        for cafa_id, e_value, bit_score in v[1]:
            query_dict.setdefault(cafa_id, {})
            query_dict[cafa_id][k] = bit_score
        gt_dict[k] = v[0]
        go_term_set = go_term_set | v[0]
    return query_dict, gt_dict, go_term_set
            
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

def parse_gaf(gaf_file, output_path):
    """
    takes a goa file, parses it and saves a the parsed version in a text file 

    Parameters
    ----------
    gaf_file : str
        path to the goa unparsed file 

    output_path : str
        path where to put the parsed goa file
    """
    uniprot_dict = dict()
    with gzip.open(gaf_file, "rb") as f:
        #unzipped_f = gzip.GzipFile(fileobj=f)
        prev_prot = ''
        for bline in f:
            line = bline.decode('UTF-8')
            if line.startswith('UniProtKB'):
                l = line.split()
                uniprot_id = l[1]
                go_id = l[4]
                namespace = l[8]
                if l[-3].startswith("taxon"):
                    taxon_id = l[-3].split(':')[1]
                elif l[-4].startswith("taxon"):
                    taxon_id = l[-4].split(':')[1]
                elif l[-5].startswith("taxon"):
                    taxon_id = l[-5].split(':')[1]
                uniprot_dict.setdefault(uniprot_id, [taxon_id, dict()])
                uniprot_dict[uniprot_id][1].setdefault(namespace, [])
                uniprot_dict[uniprot_id][1][namespace].append(go_id)
                if prev_prot == '':
                    prev_prot = uniprot_id
            if prev_prot != '' and uniprot_id != prev_prot:
                save_protein(uniprot_dict, output_path)
                uniprot_dict.pop(prev_prot, None) 
                prev_prot = uniprot_id

def save_protein(uniprot_dict, output_path):
    with open(output_path + "goa.txt", "a") as f:
        for uniprot_id, v in uniprot_dict.items():
            taxon_id = v[0]
            go_dict = v[1]
            f.write(uniprot_id + '\t' + taxon_id + '\t')
            for ns, go_list in go_dict.items():
                f.write(ns + ',')
                for idx, go_id in enumerate(go_list):
                    f.write(go_id)
                    if idx == len(go_list)-1:
                        f.write(';')
                    else:
                        f.write(',')
            f.write('\n')
    
                