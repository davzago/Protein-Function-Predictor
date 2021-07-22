#####################################
#$ -S /bin/bash
#$ -q high
#$ -cwd
#$ -pe multithread 8
##$ -e err.log
##$ -o out.log
#####################################


module load python3/3.6.12

rm -f train.sh.p{o,e}*

python3 InterPro/predict_multiple.py "InterPro/interpro_set.txt" "InterPro_output/go_dict.txt" "InterPro/models/InterPro_regression_model.joblib" "feature_matrix.txt" -output_path "InterPro/predictions"