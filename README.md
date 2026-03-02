# PspCas12f_binder
# conda environment
conda env create -f environment.yml

conda activate PspCas12f_binder
# install PyRosetta
conda install -y -c https://conda.rosettacommons.org -c conda-forge pyrosetta

conda activate CRISPRCasDesigner

cd Stragegy3-5

sh run.sh
