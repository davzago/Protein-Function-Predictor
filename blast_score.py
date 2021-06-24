
#takes the bit score dict and the ground truth dict and computes the score for each couple protein - go term
def blast_score(bit_score_dict, gt_dict, go_term_set):
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
        


    
