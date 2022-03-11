import sys
sys.path.append("../mutation_analysis")

import numpy as np
import pandas as pd
from objects.PDBData import PDBData
from objects.Mutation import Mutation
from objects.NeighborResidueSelection import NeighborResidueSelection
from features.SCONES_2021 import SCONES

# configurations
pdb_dir = "data/pdbs/"
pdbs_clean_dir = "data/pdbs_clean/"
fastas_dir = "data/fastas/"
CIF = "mmCif"

inp_file_path = "data/clean_2/PoPMuSiC_2.csv"
out_file_path = "SCONES_features_on_PoPMuSiC_2.csv"
n_rows_to_skip = 1018
n_rows_to_evalutate = 100000

# object initialization
NRS = NeighborResidueSelection()
scones = SCONES()

df = pd.read_csv(inp_file_path)

for i, row in df.iterrows():
    if i+1 <= n_rows_to_skip: continue

    # getting mutation object
    mut = Mutation().get_mutation(row)
    print(f"Row no:{i+1}->{mut}")

    cln_pdb_file = pdbs_clean_dir+mut.pdb_id+mut.chain_id+".pdb"

    neighbor_residue_ids = NRS.get_distance_based_neighbor_residue_ids(cln_pdb_file, mut.chain_id, mut.mutation_residue_id, "CB", 8.0)
    print(neighbor_residue_ids)
    for neighbor_residue_id in neighbor_residue_ids:
        features = scones.get_all_features(cln_pdb_file, mut.chain_id, mut.mutation_residue_id, neighbor_residue_id, "CB", "CB")
        print(features.shape)
        with open(out_file_path,'a') as f:
            np.savetxt(f, [features], delimiter=",")
        # break

    print()
    del mut
    if i+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break