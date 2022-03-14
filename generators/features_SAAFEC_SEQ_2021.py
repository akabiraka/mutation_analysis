import sys
sys.path.append("../mutation_analysis")

import numpy as np
import pandas as pd
from objects.Mutation import Mutation
from features.SAAFEC_SEQ_2021 import SAAFEC_SEQ
from objects.PDBData import PDBData


# configurations
# pdb_dir = "data/pdbs/"
pdbs_clean_dir = "data/pdbs_clean/"
fastas_dir = "data/fastas/"
pssms_dir = "data/pssms/"
inp_file_path = "data/clean_2/PoPMuSiC_2.csv"
out_file_path = "data/computed_features/SAAFEC_SEQ_features_on_PoPMuSiC_2.csv"
n_rows_to_skip = 0
n_rows_to_evalutate = 100000

# object initialization
saafec_seq = SAAFEC_SEQ()
pdbdata = PDBData()

df = pd.read_csv(inp_file_path)

for i, row in df.iterrows():
    if i+1 <= n_rows_to_skip: continue

    # getting mutation object
    mut = Mutation().get_mutation(row)
    print(f"Row no:{i+1}->{mut}")

    # creating necessary file paths
    pdb_file = pdbs_clean_dir+mut.pdb_id+mut.chain_id+".pdb"
    mut_pssm_file = pssms_dir+mut.pdb_id+mut.chain_id+"_"+mut.mutation_event+".pssm"
    mut_fasta_file = fastas_dir+mut.pdb_id+mut.chain_id+"_"+mut.mutation_event+".fasta"

    zero_based_mutation_site = pdbdata.get_zero_based_mutation_site(pdb_file, mut.chain_id, mut.mutation_residue_id)
    print(mut.mutation_residue_id, zero_based_mutation_site)
    features = saafec_seq.get_all_features(mut_pssm_file, mut_fasta_file, zero_based_mutation_site, 
                                           mut.wild_residue, mut.mutant_residue, window_size=7, n_neighbors=3, smooth=True)
    print(features.shape)
    with open(out_file_path,'a') as f:
        np.savetxt(f, [features], delimiter=",")
    # break

    print()
    del mut
    if i+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break