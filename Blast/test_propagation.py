import parsers


def parse_blast(blast_file):
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
            """if prev_cafa_id == '':
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
                prev_cafa_id = cafa_id"""
            if float(e_value) < 0.001:
                if cafa_id in blast_dict:
                    blast_dict[cafa_id].append([uniprot_id, float(e_value), float(bit_score)])
                else:
                    blast_dict.setdefault(cafa_id, [[uniprot_id, float(e_value), float(bit_score)]])
        #This aswell
        """rev = reverse_dict(blast_dict)
        goa_dict = parse_goa(goa_file, rev, ont_dict)
        query_dict, gt_dict, go_term_set = return_to_query(goa_dict)
        blast_score(query_dict, gt_dict, go_term_set, output_path, i)"""
    return blast_dict


goa_file = "/home/davide/Documenti/training/goa_db_2018_08_exp.dat"
ontology_file = "/home/davide/Documenti/training/gene_ontology_edit.obo.2018-08-01"
blast_file = "/home/davide/Documenti/training/blast_pred.txt"

ont, _ = parsers.parse_ontology(ontology_file)
blast_dict = parse_blast(blast_file)
rev_dict = parsers.reverse_dict(blast_dict)
i = 0
for k, v in rev_dict.items():
    if i == 10:
        break
    i += 1
    print(k,v)

uniprot_list = ["Q3TWR3", "O54950", "A0A0G2JVD6", "P80385", "G3I6J0", "A0A1S3FCV3","A0A1U7QGV1", "I3MNC0", "A0A1A6HP54", "A0A2K6SRC8"]
    
goa_dict = parsers.parse_goa(goa_file, rev_dict, ont)
for k in uniprot_list:
    if k in goa_dict:
        print(k, goa_dict[k])

"""go_set = {'GO:0051170', 'GO:0005829', 'GO:0004679', 'GO:0005634', 'GO:0005654'}
print(go_set)
prop_set = parsers.propagate_terms({'GO:0051170', 'GO:0005829', 'GO:0004679', 'GO:0005634', 'GO:0005654'}, ont)
print(prop_set)"""
