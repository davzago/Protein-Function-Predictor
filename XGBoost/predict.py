import xgboost as xgb
import argparse
from utils import *

parser = argparse.ArgumentParser(description='XGBoost predictor')
parser.add_argument('Naive_predictions', help='Path to the folder containing the predictions from the Naive component')
parser.add_argument('Blast_predictions', help='Path to the folder containing the predictions from the Blast component')
parser.add_argument('InterPro_predictions', help="Path to the folder containing the predictions from the Interpro component")
parser.add_argument('model', help="model used to make the predictions")
parser.add_argument('k', help="maximum number of proposed go term from each component")
parser.add_argument('-output_path', help='Path to the folder where the predction file will be put', default="output")
args = parser.parse_args()

naive_folder = args.Naive_predictions
blast_folder = args.Blast_predictions
interpro_folder = args.InterPro_predictions

model = args.model
k = args.k
output_path = args.output_path

rank_model = xgb.Booster({'nthread': 8})
rank_model.load_model(model)

naive_pred = parse_component_prediction(naive_folder)
blast_pred = parse_component_prediction(blast_folder)
interpro_pred = parse_component_prediction(interpro_folder)

pred_dict = combine_dictionaries(naive_pred, blast_pred, interpro_pred, k)

df, assoc = build_dataset(pred_dict)

X = df[['Naive', 'Blast', 'InterPro']].to_numpy().astype('float')

groups = df.groupby('Protein id').size().to_frame('size')['size'].to_numpy()

pred = predict_groups(rank_model, X, groups)

prot_ids = df[['Protein id']].to_numpy()
go_terms = df[['Go term id']].to_numpy()

save_predictions(pred, prot_ids, go_terms)

