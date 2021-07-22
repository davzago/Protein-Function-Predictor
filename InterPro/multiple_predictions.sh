#####################################
#$ -S /bin/bash
#$ -q high
#$ -cwd
##$ -e err.log
##$ -o out.log
#####################################

line=$(sed "${SGE_TASK_ID}q;d" $1)

module load python3

rm -f multiple_predictions.sh.p*

python3 InterPro/predict.py "interpro_set.txt" "InterPro_output/go_dict.txt" "InterPro/models/InterPro_regression_model.joblib" "$line" -output_path "InterPro_predictions"
