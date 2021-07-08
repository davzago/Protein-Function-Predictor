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