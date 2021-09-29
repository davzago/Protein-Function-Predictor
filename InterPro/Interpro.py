from parsers import parse_dict, parse_dict_to_array, parse_features
from dataset_utils import *
from dataset_utils import parse_interpro_set, reverse_dict, save_prediction, set_ip_indices
from joblib import load


def interpro_predict(model, interpro_set, terms_dict, features):
    """
    Takese a set of proteins and assigns to them go terms based on their interpro labels

    Parameters
    ----------
    model : 
        a multiple output classification model
    
    interpro_set : str
        file containing the interpro labels found during the training of the model

    terms_dict : str
        file that contains the mapping between the go terms and their position in the output

    features : str
        tsv file containing the proteins and their interpro label

    Returns
    -------
    pred_dict : dict()
        A dictionary containing the predictions done by the model {prot_id : [[go_term, score]]}
    """

    clf = load(model)
    ip_set = parse_interpro_set(interpro_set)
    #go_dict = parse_dict(terms_dict)
    #go_array = reverse_dict(go_dict)
    go_array = parse_dict_to_array(terms_dict)

    X, ip_dict, prot_dict = parse_features(features, ip_set)

    pred = clf.predict_proba(X)

    prot_array = reverse_dict(prot_dict)

    Y = reshape_prediction(pred)

    pred_dict = prediction_to_dict(Y, prot_array, go_array)

    return pred_dict