import os
import argparse

parser = argparse.ArgumentParser(description='Blast predictor')
parser.add_argument('uniprot_file', help='File containing the uniprot ids you want to find in the fasta')
parser.add_argument('fasta_file', help="File that contains all the matches between cafa id and uniprot id")
parser.add_argument('file_name', help="name of the file")
args = parser.parse_args()


uniprot_file = args.uniprot_file
fasta_file = args.fasta_file
filename = args.file_name

def get_set(uniprot_list):
    uniprot_set = set()
    with open(uniprot_list, 'r') as f:
        for line in f:
            uniprot_id = line.split()[0]
            uniprot_set.add(uniprot_id)
    return uniprot_set

def parse_fasta(uniprot_set, fasta_file, filename, n_max):
    n = 0
    n_file = 1
    uniprot_id = ''
    with open(fasta_file, 'r') as fasta:
        sequence = ''
        for line in fasta:
            if line.startswith('>'):
                if sequence != '' and uniprot_id in uniprot_set:
                    n += 1
                    with open(filename + str(n_file) + ".fasta", 'a') as f:
                        f.write(">" + uniprot_id + '\n' + sequence)
                sequence = ''
                uniprot_id = line.split('|')[1] 
                if n >= n_max:
                    n = 0
                    n_file += 1
            elif line.isupper():
                sequence += line.strip()


uniprot_set = get_set(uniprot_file)
parse_fasta(uniprot_set, fasta_file, filename , 10000)