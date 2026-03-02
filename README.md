# PspCas12f_binder
# conda environment & install PyRosetta
```bash
conda env create -f environment.yml
conda activate PspCas12f_binder
conda install -y -c https://conda.rosettacommons.org -c conda-forge pyrosetta
```
# Strategy1:Homology-based design
```bash
cd Stragegy1
python Blast_feature.py
```
# Strategy3-5:Structure-informed rational design
```bash
cd Strategy3-5
name=fold_pspcas12f1_a16_model_0
python3.10 script/find_bind.py pdb_af
python3.10 script/find_bind_dimer.py pdb_af A B 5
python3.10 script/pdb_mutation_all.py mutlist.txt pdb_af change_result/$name
python3.10 script/find_bind_muti.py change_result/$name
python3.10 script/bind_compare.py
cp pdb_af/$name.pdb change_result/$name/
python3.10 script/rosetta_all.py change_result/$name
python3.10 script/bonds_select.py change_result/${name}_bonds $name
```
