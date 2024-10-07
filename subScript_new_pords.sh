#!/bin/bash

#SBATCH --nodes=1

#SBATCH --time=05:30:00

#SBATCH --mem=300gb
 
#SBATCH --partition=sla-prio

#SBATCH --account=cdm8


# Get started


echo "Here we go"
echo "--------------"
#module load cplex
#module load anaconda3
#sconda init bash
#conda activate chemicalrxn
#pip install rdkit-pypi
#cd $SBATCH_O_WORKDIR
#cd /storage/group/cdm8/default/mohit/chemicalrxn/Retrosynthesis_planning/py_data/new_data/
#cd ./../batch_code/

echo "All the code starts here"
echo "------------------------"

#python milp4_Hydroxy_butyrate_glu.py
#python milp4_3HB_gly.py
#python milp4_3HB_pyr.py
#python milp4_Limonene_R_pyr.py
#python milp4_Limonene_R_pXy.py
#python milp4_Tereph_pXy.py
#python milp4_FDCA_glu.py
#python milp4_CHEMBL502135_glu.py
#python milp4_CHEMBL282575_pXy.py
#python milp4_CHEMBL590540_pyr.py
python minchembio.py







# Finish up 


