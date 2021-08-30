import argparse

def get_blast_pids(blast_file):
    """
    parses the blast file and returns a list of protein ids

    Parameters
    ----------
    blast_file: str 
        path to the file containing the goa db 

    Returns
    -------
    cafa : set
        set containing protein's cafa ids
    """
    cafa = set()
    with open(blast_file, "r") as f:
        for line in f:
            div = line.split()
            cafa.add(div[0])
    return cafa

def ref_to_benchmark(ref_file):
    """
    Takese a reference file and returns the list of cafa ids in it

    Parameters
    ----------
    ref_file : str
        path to the reference file
    """

    benchmark = dict()
    with open(ref_file, 'r') as f:
        for line in f:
            cafa_id, go_id, namespace, uniprot_id = line.split()
            benchmark.setdefault(namespace, set())
            benchmark[namespace].add(cafa_id)
    return benchmark

def save_benchmark(benchmark_dict, name, output_folder):
    for ns, set in benchmark_dict.items():
        if ns == 'F':
            with open(output_folder + '/' + name + '_mfo.txt', 'w') as f:
                for cafa_id in set:
                    f.write(cafa_id + '\n')
        elif ns == 'C':
            with open(output_folder + '/' + name + '_cco.txt', 'w') as f:
                for cafa_id in set:
                    f.write(cafa_id + '\n')
        elif ns == 'P':
            with open(output_folder + '/' + name + '_bpo.txt', 'w') as f:
                for cafa_id in set:
                        f.write(cafa_id + '\n')        
             

parser = argparse.ArgumentParser(description='Outputs a benchmark file starting from a reference file')
parser.add_argument('reference_file', help='File containing the reference')
parser.add_argument('blast_file', help='File containing the Blast results')
parser.add_argument('name', help='name of the output file')
parser.add_argument('-output_path', help='Path to the folder where the benchmark file will be put', default="output")
args = parser.parse_args()

ref = args.reference_file
blast = args.blast_file
name = args.name
output = args.output_path


benchmark = ref_to_benchmark(ref)
blast_set = get_blast_pids(blast)
for ns, set in benchmark.items():
    benchmark[ns] = set.intersection(blast_set) 
save_benchmark(benchmark, name, output)