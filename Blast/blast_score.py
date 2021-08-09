
# takes the bit score dict and the ground truth dict and computes the score for each couple protein - go term
# score_dict = {cafa_id: {go_term: score}}
def blast_score(bit_score_dict, gt_dict, go_term_set, output_path, i):
    """
    takes the bit score dict and the ground truth dict and computes the score for each pair of proteins

    Parameters
    ----------
    bit_score_dict: dict
        {query_id: {uniprot_id: bit_score}}

    gt_dict : dict
        {uniprot_id: go_term_set}
    
    go_term_set : set
        set containing all the go terms found in the goa database

    output_path : str
        path to the output folder

    i : int 
        parameter used to name the various prediction files
    """
    score_dict = dict()
    for k, v in bit_score_dict.items():
        score_dict.setdefault(k, dict())
        for term in go_term_set:
            n = 0
            d = 0
            for uniprot_id, bit_score in v.items():
                d += bit_score
                if term in gt_dict[uniprot_id]:
                    n += bit_score
            score_dict[k][term] = n/d 
    prediction_to_text(score_dict, "prediction"+str(i), output_path)
    score_dict = dict()
    #return score_dict

def blast_score_no_out(bit_score_dict, gt_dict, go_term_set):
    """
    takes the bit score dict and the ground truth dict computes the score for each pair of proteins, returns  

    Parameters
    ----------
    bit_score_dict: dict
        {query_id: {uniprot_id: bit_score}}

    gt_dict : dict
        {uniprot_id: go_term_set}
    
    go_term_set : set
        set containing all the go terms found in the goa database

    Returns
    -------
    score_dict : dict
        dict : {query_id : {go_term : blast_score}}
    """
    score_dict = dict()
    for k, v in bit_score_dict.items():
        score_dict.setdefault(k, dict())
        for term in go_term_set:
            n = 0
            d = 0
            for uniprot_id, bit_score in v.items():
                d += bit_score
                if term in gt_dict[uniprot_id]:
                    n += bit_score
            score_dict[k][term] = n/d 
    return score_dict
        
def prediction_to_text(score_dict, file_name, output_path='', n_go_terms=-1):
    """
    takes a dictionary containing the score for each go term corresponding to the prot_id
    and writes the prediction files for the evaluation

    Parameters
    ----------
    score_dict : dict
        {query_id : {go_id: score}}

    file_name : str
        name of the prediction file where the data will be written

    output_path : str
        path to the output folder

    n_go_terms : int
        maximum number of go terms to associate with the protein    
    """
    f = open(output_path +'/'+ file_name + ".txt", "w")
    for cafa_id, go_dict in score_dict.items():
        i = 0
        for go_id, score in sorted(go_dict.items(), key=lambda x: x[1], reverse=True):
            if score > 0:
                f.write(cafa_id + "\t" + go_id + "\t" + str(score) + "\n")
            if i == n_go_terms:
                break
            i += 1
    f.close()


def save_predictions(score_dict, file_name ='predictions', output_path='', n_go_terms=-1):
    """
    takes a dictionary containing the score for each go term corresponding to the prot_id
    and writes the prediction files for the evaluation

    Parameters
    ----------
    score_dict : dict
        {query_id : {go_id: score}}

    file_name : str
        name of the prediction file where the data will be written

    output_path : str
        path to the output folder

    n_go_terms : int
        maximum number of go terms to associate with the protein    
    """
    count = 0
    j = 1
    f = open(output_path +'/'+ file_name + str(j) + ".txt", "w")
    for cafa_id, go_dict in score_dict.items():
        if count >= 500:
            f.close()
            j += 1
            f = open(output_path +'/'+ file_name + str(j) + ".txt", "w")
            count = 0
        i = 0
        for go_id, score in sorted(go_dict.items(), key=lambda x: x[1], reverse=True):
            if score > 0:
                f.write(cafa_id + "\t" + go_id + "\t" + str(score) + "\n")
            if i == n_go_terms:
                break
            i += 1
        count += 1
    f.close()

