#####################################
#$ -S /bin/bash
#$ -q high
#$ -cwd
##$ -e err.log
##$ -o out.log
#####################################

line=$(sed "${SGE_TASK_ID}q;d" $1)

module load python3

rm -f predictions/*.txt

python3 InterPro/predict.py "InterPro/interpro_set.txt" "InterPro/go_dict.txt" "InterPro_output/model.joblib" $line -output_path "predictions"
