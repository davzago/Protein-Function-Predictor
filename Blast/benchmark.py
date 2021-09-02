import argparse
import os
from utils import *
      
             
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

if not os.path.isdir(output):
    os.mkdir(output)

benchmark = ref_to_benchmark(ref)
blast_set = get_blast_pids(blast)
for ns, set in benchmark.items():
    benchmark[ns] = set.intersection(blast_set) 
save_benchmark(benchmark, name, output)