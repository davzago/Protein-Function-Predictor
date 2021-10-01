import os

uniprot_list = "uniprot_list.txt"
results_folder = "blast_results"

result_list = os.listdir(results_folder)

uniprot_set = set()
with open(uniprot_list, 'r') as f:
    for line in f:
        uni_id = line.strip()
        uniprot_set.add(uni_id)

for result in result_list:
    if result.endswith("fasta"):
        with open(results_folder + '/' + result, 'r') as f:
            for line in f:
                uni_id = line.split()[0]
                if uni_id in uniprot_set:
                    uniprot_set.remove(uni_id)

with open("updated_uniprot_list.txt", 'w') as f:
    for member in uniprot_set:
        f.write(member + '\n')
