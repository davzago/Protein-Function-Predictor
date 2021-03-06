import os
import numpy as np
import pandas as pd
import xgboost as xgb
import math

def parse_prediction(prediction_file, prediction_dictionary, n_component, k):
    """
    Takes a prediction file and returns a dictionary containing the prot id as key and a dict of scores for each go term as value

    Parameters
    ----------
    prediction_file : str
        path to the prediction file 

    prediction_dictionary : dict
        build in progress dictionary that will contain all the predictions: {cafa_id: {go_term: [scores]}}

    n_component : int
        the index where the score of the component should be placed
    
    k : int
        maximum number of go terms per protein
    Returns
    -------
    prediction_dictionary : dict
        same dict taken as input augmented with the predictions found in the prediction file 
    """
    component_dict = dict()
    with open(prediction_file, 'r') as f:
        for line in f:
            cafa_id, go_term, score = line.split()
            component_dict.setdefault(cafa_id, [])
            component_dict[cafa_id].append([go_term, score])

    for cafa_id, go_list in component_dict.items():
        go_list.sort(key=lambda x: float(x[1]))
        for go_term, score in go_list[::-1][:k]:
            prediction_dictionary.setdefault(cafa_id, dict())
            prediction_dictionary[cafa_id].setdefault(go_term, [[0,0,0],0])
            prediction_dictionary[cafa_id][go_term][0][n_component] = score

    return prediction_dictionary

def build_prediction_dict(prediction_folder_list, k):
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
    for i, folder in enumerate(prediction_folder_list):
        for pred_file in os.listdir(folder):
             pred_dict = parse_prediction(folder +'/'+ pred_file, pred_dict, i, k)
    return pred_dict

def parse_component_prediction(prediction_folder):
    """
    Takes the folder containing all the prediction of a component and build a dictionary with them

    Parameters
    ----------
    prediction_folder : str
        path to the folder containing the predictions

    Returns
    -------
    prediction_dict : dict()
        {prot_id : {go_id : score}}
    """
    pred_dict = dict()
    for pred_file in os.listdir(prediction_folder):
        if pred_file.startswith("pred"):
            with open(prediction_folder + '/' + pred_file, 'r') as f:
                for line in f:
                    cafa_id, go_term, score = line.split()
                    pred_dict.setdefault(cafa_id, dict())
                    pred_dict[cafa_id][go_term] = score
    return pred_dict

def combine_dictionaries(pred1, pred2, pred3, k):
    """
    takes the predictions dictionaries of the 3 components, the number of proposed terms and builds the dataframe

    Parameters
    ----------
    pred1: dict()
        {prot_id : {go_id : score}}

    pred1: dict()
        {prot_id : {go_id : score}}
    
    pred1: dict()
        {prot_id : {go_id : score}}
    
    k : int
        maximum number of proposed terms for each protein

    Returns
    -------
    combined_dict : dict()
        {protein_id: {go_id: [[score1, score2, score3], gt]}}
    """
    final_dict = dict()
    for prot_id , go_dict in pred1.items():
        go_list = []
        for go_id, score in go_dict.items():
            go_list.append([go_id, score])
        go_list.sort(key=lambda x: float(x[1]))
        go_list = go_list[::-1][:k]
        final_dict.setdefault(prot_id, dict())
        for go, s in go_list:
            final_dict[prot_id][go] = [[float(s), 0, 0], 0]

    for prot_id, go_dict in pred2.items():
        go_list = []
        for go_id, score in go_dict.items():
            go_list.append([go_id, score])
            if prot_id in final_dict and go_id in final_dict[prot_id]:
                final_dict[prot_id][go_id][0][1] = score
        go_list.sort(key=lambda x: float(x[1]))
        go_list = go_list[::-1][:k]
        final_dict.setdefault(prot_id, dict())
        for go, s in go_list:
            if go not in final_dict[prot_id]:
                final_dict[prot_id][go] = [[0, float(s), 0], 0]

    for prot_id, go_dict in pred3.items():
        go_list = []
        for go_id, score in go_dict.items():
            go_list.append([go_id, score])
            if prot_id in final_dict and go_id in final_dict[prot_id]:
                final_dict[prot_id][go_id][0][2] = score
        go_list.sort(key=lambda x: float(x[1]))
        go_list = go_list[::-1][:k]
        final_dict.setdefault(prot_id, dict())
        for go, s in go_list:
            if go not in final_dict[prot_id]:
                final_dict[prot_id][go] = [[0, 0, float(s)], 0]

    return final_dict
                
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
    df : pandas DataFrame
        Dataframe containing dataset and labels columns=['Protein id', 'group id', 'Go term id', 'Naive', 'Blast', 'InterPro', 'Label']

    assoc : dict
        dictionary that stores the which group corresponds to each protein
    """
    assoc = dict()
    next_free = 0
    data = []
    

    for prot, go_dict in pred_dict.items():
        #fill the association list
        if prot not in assoc:
            assoc[prot] = next_free
            next_free += 1
        #group_size = 0
        for go_term, [[score1, score2, score3], label] in go_dict.items():
            data.append([prot, assoc[prot], go_term, score1, score2, score3, label])
            #group_size += 1
    df = pd.DataFrame(data, columns=['Protein id', 'Group', 'Go term id', 'Naive', 'Blast', 'InterPro', 'Label'])
    return df, assoc

def add_ground_truth_from_goa(goa, pred_dict): 
    """
    takes a prediction dict and adds the ground truth contained in the goa file

    Parameters
    ----------
    goa : str
        goa file containing the ground truth

    pred_dict : dict
        the dictionary containing the predictions {prot_id: {go_id: [[score1,score2,score3],gt]}}
    """

    with open(goa) as f:
        for line in f:
            uniprot_id, _, terms = line.split()
            for go_list in terms.split(';'):
                l = go_list.split(',')
                for go_term in l:
                    if uniprot_id in pred_dict and go_term in pred_dict[uniprot_id]:
                        pred_dict[uniprot_id][go_term][1] = 1
    return pred_dict

def add_ground_truth(ref, pred_dict): 
    # this should be changend for the final version
    """
    takes a prediction dict and adds the ground truth contained in the ref file

    Parameters
    ----------
    ref : str
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
            if cafa_id in pred_dict and go_term in pred_dict[cafa_id]:
                pred_dict[cafa_id][go_term][1] = 1
    return pred_dict
    
def predict(model, df):
        return model.predict(df.iloc[:, ~df.columns.isin(['Group'])])

def predict_groups(model, data, groups):
    """
    takes a model, the features and the groups information and computes the prediction and normalizes them in the [0,1] interval

    Parameters
    ----------
    model : trained model
        the ranking model that will be use to make predictions

    data : numpy matrix
        matrix that has proteins as rows and components score as columns

    groups : numpy array
        list of integers representing the size of subsequent group of feature (there is 1 group for each protein)

    Returns
    -------
    pred : numpy array
        list of prediction scores computed using the model
    """
    pred = []
    index = 0
    for g in groups:
        dX = xgb.DMatrix(data[index:index+g,:])
        prediction = model.predict(dX)
        prediction = (prediction - min(prediction)) / (max(prediction) - min(prediction))
        pred = [*pred, *prediction]
        index += g
    #pred = (pred - min(pred)) / (max(pred) - min(pred))
    return np.array(pred)

def scores_mean(list_of_scores):
    mean = sum(list_of_scores)/len(list_of_scores)
    return mean

def save_predictions(preds, prot_ids, go_terms, groups, n_prot, output_path):
    """
    saves the predictions made by the xgboost component in a tsv file 

    Parameters
    ----------
    preds : numpy array
        array containing the predictions scores

    prot_ids : numpy array 
        array containing the protein's id to match each score in preds

    go_terms : numpy array
        array containing the go terms to match each score in preds 

    groups . numpy array
        array containing how many elements are in each group (same protein id)

    n_prot : int
        maxium number of proteins per prediction file

    output_path : str
        path where the prediction file will be put
    """
    go_len = len(go_terms)
    prot_len = len(prot_ids)
    preds_len = len(preds)
    n = 1
    div = math.ceil(len(groups) / n_prot) 
    j = 0
    if prot_len == preds_len and go_len == prot_len:
        for i in range(0,div):
            start = i * n_prot
            end = start + n_prot 
            g = sum(groups[start:end])
            e = j + g
            m = np.array([prot_ids[j:e], go_terms[j:e], preds[j:e]]).T
            test = set(prot_ids[j:e])
            np.savetxt(output_path + '/' + "predictions"+ str(n) +".txt", m, delimiter='\t', fmt='%s')
            j = e
            n += 1
    else:
        raise Exception("the provided arrays do not have the same length")

def switch_prot_ids(mapping_file, prot_ids):
    """
    takes a list of uniprot ids and converts them into a list of cafa ids

    Parameters
    ----------
    mapping_file: str
        path to the file that contains the mapping between uniprot ids and cafa ids

    prot_ids: numpy array
        list of uniprot ids

    Returns
    -------
    prot_list: numpy array
        list of cafa ids (will be used to save the predictions)
    """       
    id_map = dict()
    prot_list = []
    with open(mapping_file, 'r')as f:
        for line in f:
            cafa_id, uniprot_id , _ = line.split()
            id_map.setdefault(uniprot_id, cafa_id)
    
    for i in range(0,len(prot_ids)):
        if prot_ids[i] in id_map:
            prot_list.append(id_map[prot_ids[i]])
    
    return np.array(prot_list)
