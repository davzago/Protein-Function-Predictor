
# takes the bit score dict and the ground truth dict and computes the score for each couple protein - go term
# score_dict = {cafa_id: {go_term: score}}
def blast_score(bit_score_dict, gt_dict, go_term_set):
    score_dict = dict()
    counter = dict
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
        prediction_to_text(score_dict, k, "predictions/Fake")
        score_dict = dict()
    #return score_dict
        
def prediction_to_text(score_dict, prot_id, output_path='', n_go_terms=-1):
    f = open(output_path +'/'+ prot_id + ".txt", "a")
    for cafa_id, go_dict in score_dict.items():
        i = 0
        for go_id, score in sorted(go_dict.items(), key=lambda x: x[1], reverse=True):
            f.write(cafa_id + "\t" + go_id + "\t" + str(score) + "\n")
            if i == n_go_terms:
                break
            i += 1
    f.close()

