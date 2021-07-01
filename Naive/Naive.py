
# Parses the goa file and retrurns a dictionary dict{}
def parse_goa(goa_file):
    """
    Parses a goa file and returns its data as a dictionary

    Parameters
    ----------
    goa_file : str
        the file location of the goa_db

    Returns
    -------
    dictionary
        {uniprot_id: set(go_id)}
    """
    goa_dict = dict()
    with open(goa_file) as f:
        for line in f:
            uniprot_id, _, terms = line.split()
            for go_list in terms.split(';'):
                l = go_list.split(',')
                for i in range(1,len(l)):
                    if uniprot_id in rev_dict:
                        goa_dict.setdefault(uniprot_id, set())
                        goa_dict[uniprot_id].add(l[i])
    return goa_dict


def Naive_calculator(goa_dict):
    """
    Computes the frequency of go terms in a dataset

    Parameters
    ----------
    goa_dict : dictionary
        {uniprot_id: set(go_id)}

    Returns
    -------
    dictionary
        {go_id: frequency}
    """
    
