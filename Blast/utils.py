def get_reference_pids(ref_file):
    """
    parses the ref file and returns a list of protein ids

    Parameters
    ----------
    ref_file: str 
        path to the file containing the reference 

    Returns
    -------
    uniprot : set
        set containing protein ids

    cafa : set
        set containing protein ids
    """
    uniprot = set() 
    cafa = set()
    with open(ref_file, "r") as f:
        for line in f:
            div = line.split()
            uniprot.add(div[-1])
            cafa.add(div[0])
    return uniprot, cafa

def get_goa_pids(goa_file):
    """
    parses the ref file and returns a list of protein ids

    Parameters
    ----------
    goa_file: str 
        path to the file containing the goa db 

    Returns
    -------
    uniprot : set
        set containing protein ids
    """
    uniprot = set()
    with open(goa_file, "r") as f:
        for line in f:
            div = line.split()
            uniprot.add(div[0])
    return uniprot

def get_blast_pids(blast_file):
    """
    parses the ref file and returns a list of protein ids

    Parameters
    ----------
    goa_file: str 
        path to the file containing the goa db 

    Returns
    -------
    uniprot : set
        set containing protein ids
    """
    cafa = set()
    with open(blast_file, "r") as f:
        for line in f:
            div = line.split()
            cafa.add(div[0])
    return cafa

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


