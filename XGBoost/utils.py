import os

def parse_prediction(prediction_file, prediction_dictionary):
    """
    Takes a prediction file and returns a dictionary containing the prot id as key and a dict of scores for each go term as value

    Parameters
    ----------
    prediction_file : str
        path to the prediction file 

    prediction_dictionary : dict
        build in progress dictionary that will contain all the predictions: {cafa_id: {go_term: [scores]}}
    
    Returns
    -------
    prediction_dictionary : dict
        same dict taken as input augmented with the predictions found in the prediction file 
    """
    
    with open(prediction_file, 'r') as f:
        for line in f:
            cafa_id, go_term, score = line.split()
            prediction_dictionary.setdefault(cafa_id, dict())
            prediction_dictionary[cafa_id].setdefault(go_term, [[],-1])
            prediction_dictionary[cafa_id][go_term].append(float(score))
    return prediction_dictionary


def build_prediction_dict(prediction_folder_list):
    """
    Takes the folders containing the predictions and outputs a dictionary that groups all the predictions

    Parameters
    ----------
    prediction_folder_list : list
        list of prediction folders one for each method

    Returns
    -------
    prediction_dictionary : dict
        dictionary containing all the predictions
    """
    pred_dict = dict()
    for folder in prediction_folder_list:
        for pred_file in os.listdir(folder):
             pred_dict = parse_prediction(pred_file, pred_dict)
    return pred_dict

def build_dataset(pred_dict):
    """
    takes the predictions of the feature extraction components of the predictor as a dictionary and returns them
    in the right format to build the dataset

    Parameters
    ----------
    pred_dict : dict
        {prot_id: {go_id: [[score1,score2,score3],gt]}}

    Returns
    -------
    feature_matrix : np_array
        matrix containing the scores for each go term

    labels : np_array
        array containing the ground truth

    group : np_array 
        array containing the group id for each line in the feature matrix

    assoc : dict
        dictionary that stores the which group corresponds to each protein
    """
    assoc = dict()
    next_free = 0

    for prot, scores in pred_dict.items():
        if prot not in assoc:
            assoc[prot] = next_free
            next_free += 1
        
#def add_ground_truth(goa, pred_dict):
    """
    takes a prediction dict and adds the ground truth contained in the goa file

    Parameters
    ----------
    goa : str
        file containing the ground truth

    pred_dict : dict
        the dictionary containing the predictions {prot_id: {go_id: [[score1,score2,score3],gt]}}
    """

def add_ground_truth(ref, pred_dict):
    """
    takes a prediction dict and adds the ground truth contained in the ref file

    Parameters
    ----------
    goa : str
        file containing the ground truth

    pred_dict : dict
        the dictionary containing the predictions {prot_id: {go_id: [[score1,score2,score3],gt]}}

    Returns
    -------
    pred_dict : dict
        pred_dict taken in imput with the addition of the ground truth
    """

    with open(ref, 'r') as f:
        for line in f:
            cafa_id, go_term, namespace, uniprot_id = line.split('\t')
            if cafa_id in pred_dict:
                pred_dict[]

        
