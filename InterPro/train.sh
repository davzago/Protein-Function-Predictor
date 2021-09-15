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
python3 InterPro/train.py "interpro_set.txt" "training/gene_ontology_edit.obo.2018-08-01" "goa_dataset.txt" -output_path "InterPro_output" -model_name "full_goa_2018_model"
