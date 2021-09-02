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
            prediction_dictionary[cafa_id].setdefault(go_term, [])
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