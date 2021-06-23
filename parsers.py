import pandas as pd

# could divide the list of ancestors in namespaces
def parse_ontology(obo_file_path):
    go_terms = dict()
    alts = dict()
    temp_id = ''
    temp_namespace = '' 
    temp_name = ''
    temp_def = ''
    temp_alt_id = []
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
                elif k == "name":
                    temp_name = v
                elif k == "namespace" and v != 'external':
                    temp_namespace = v
                elif k == "def":
                    temp_def = v
                elif k == "is_a" and (v.startswith("GO:") or v.startswith("HP:") or v.startswith("DO:")):
                    # add (temp_id,v) tuple to the relation list
                    s = v.split('!')[0].strip()
                    rel_list.append(s)

        go_terms[temp_id] = rel_list
        alts[temp_id] = temp_alt_id

    return go_terms, alts

def parse_blast(blast_file, prot_id):
    data = pd.read_csv(blast_file, sep='\t', names=['cafa_id', 'uniprot_id', 'identity', 'alignment_leght',
     'mismatch', 'gap_opens', 'query_start', 'query_end', 'sub_start', 'sub_end', 'e_value', 'bit_score'])
    
    blast_dict = dict()

    fil = data["cafa_id"] == prot_id
    for row in data[fil].loc[:, ['uniprot_id', 'e_value', 'bit_score']].itertuples():
        blast_dict[row[1]] = {'e_value': row[2], 'bit_score': row[3]}
        

    return blast_dict

parse_blast("/home/davide/Documenti/training/blast_pred.txt", "T100900000100")

def parse_goa(goa_file):