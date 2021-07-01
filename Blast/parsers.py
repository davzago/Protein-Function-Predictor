import pandas as pd
import os
from blast_score import *
import time

# could divide the list of ancestors in namespaces
def parse_ontology(obo_file_path):
    go_terms = dict()
    alts = dict()
    temp_id = ''
    temp_namespace = '' 
    temp_name = ''
    temp_def = ''
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

# takes a reference file and splits it in a file for each namespace
# (the uniport ids are also removed because they aren't useful to the evaluation)
def parse_reference(reference_file):
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


#takes a blast output and returns the corresponding python dictionary
def parse_blast(blast_file):
    """data = pd.read_csv(blast_file, sep='\t', names=['cafa_id', 'uniprot_id', 'identity', 'alignment_leght',
     'mismatch', 'gap_opens', 'query_start', 'query_end', 'sub_start', 'sub_end', 'e_value', 'bit_score'])
    
    blast_dict = dict()

    for row in data.loc[:, ['cafa_id', 'uniprot_id', 'e_value', 'bit_score']].itertuples():
        blast_dict.setdefault(row[1], [])
        if row[3] < 0.001:
            blast_dict[row[1]].append([row[2], row[3], row[4]])
    return blast_dict"""
    blast_dict = dict()
    with open(blast_file) as f:
        prev_cafa_id = ''
        count = 0
        for line in f:
            cafa_id, uniprot_id, _, _, _, _, _, _, _, _, e_value, bit_score = line.split()
            ### When having memory problems uncomment ths code and run only this function (will write the prediction file)
            if prev_cafa_id == '':
                prev_cafa_id = cafa_id
            if prev_cafa_id != cafa_id:
                count += 1
                if count > 10000:
                    rev = reverse_dict(blast_dict)
                    goa_dict = parse_goa("fake_goa.txt", rev, ont_dict)
                    query_dict, gt_dict, go_term_set = return_to_query(goa_dict)
                    blast_score(query_dict, gt_dict, go_term_set)
                    count = 0
                prev_cafa_id = cafa_id
            if float(e_value) < 0.001:
                if cafa_id in blast_dict:
                    blast_dict[cafa_id].append([uniprot_id, float(e_value), float(bit_score)])
                else:
                    blast_dict.setdefault(cafa_id, [[uniprot_id, float(e_value), float(bit_score)]])
        #This aswell
        rev = reverse_dict(blast_dict)
        goa_dict = parse_goa("fake_goa.txt", rev, ont_dict)
        query_dict, gt_dict, go_term_set = return_to_query(goa_dict)
        blast_score(query_dict, gt_dict, go_term_set)
    return blast_dict

# takes a blast dict b_dict={cafa_id: [uniprot_id, e_value, bit_score]} and returns a dictionary
# using uniprot_id as key uni_dict = {uniprot_id: (go_terms_set,[cafa_id, e_value, bit_score])}
def reverse_dict(blast_dict):
    rev_dict = dict()
    for cafa_id in list(blast_dict):
        for uniprot_id, e_value, bit_score in blast_dict[cafa_id]:
            rev_dict.setdefault(uniprot_id, [set(),[]])
            rev_dict[uniprot_id][1].append([cafa_id, e_value, bit_score])
        del blast_dict[cafa_id]
    return rev_dict

#Parses a GOA file and adds to the uniprot_dict(obtained in reverse_dict) a set of associated go-terms
def parse_goa(goa_file, rev_dict, ont_dict):
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
    queue = list(go_term_set)
    while queue:
        term = queue.pop(0)
        if term in ont:
            parents = ont[term]
            for p in parents:
                go_term_set.add(p)
                queue.append(p)
    return go_term_set

start = time.time()         
ont_dict, _ = parse_ontology("/home/davide/Documenti/training/gene_ontology_edit.obo.2018-08-01")
blast_dict = parse_blast("fake_blast.txt")
"""print("blast parsed")
rev = reverse_dict(blast_dict)
blast_dict.clear()
print("dictionary reversed")
goa_dict = parse_goa("/home/davide/Documenti/training/goa_db_2018_08_exp.dat", rev, ont_dict)
rev.clear()
print("goa parsed")
query_dict, gt_dict, go_term_set = return_to_query(goa_dict)
goa_dict.clear()
print("back to query")
score_dict = blast_score(query_dict, gt_dict, go_term_set)
print("score calculation")
prediction_to_text(score_dict, 100)"""
#parse_reference("/home/davide/Documenti/training/reference_new.txt")


#print(blast_dict['T100900012018'])

"""i=0
for k, v in sorted(score_dict['T100900012018'].items(), key=lambda p:p[1], reverse=True):
    print(k, v)
    if i == 20:
        break
    i += 1"""

print("time elapsed:",time.time()-start)

