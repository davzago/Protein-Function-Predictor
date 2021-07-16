#####################################
#$ -S /bin/bash
#$ -q high
#$ -cwd
##$ -e err.log
##$ -o out.log
#####################################


module load python3/3.6.12



rm -f train.sh.{o,e}*
python3 InterPro/train.py "interpro_set.txt" "training/gene_ontology_edit.obo.2018-08-01" "dataset.txt" -output_path "InterPro_output"
