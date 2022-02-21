import sys
sys.path.append("../mutation_analysis")

import csv
import pandas as pd
import numpy as np
import torch
from objects.PDBData import PDBData
from objects.Selector import ChainAndAminoAcidSelect
import generators.utils as Utils
from objects.Mutation import Mutation

# configurations
pdb_dir = "data/pdbs/"
pdbs_clean_dir = "data/pdbs_clean/"
fastas_dir = "data/fastas/"
CIF = "mmCif"

input_file_path = "data/clean_1/PoPMuSiC_2.csv"
n_rows_to_skip = 0
n_rows_to_evalutate = 10#0000

# object initialization
pdbdata = PDBData(pdb_dir=pdb_dir)

df = pd.read_csv(input_file_path)

for i, row in df.iterrows():
    if i+1 <= n_rows_to_skip: continue

    # getting mutation object
    mut = Mutation().get_mutation(row)
    # print(mut)
  
    # creating necessary file paths
    cln_pdb_file = pdbs_clean_dir+mut.pdb_id+mut.chain_id+".pdb"
    wild_fasta_file = fastas_dir+mut.pdb_id+mut.chain_id+".fasta"
    mutant_fasta_file = fastas_dir+mut.pdb_id+mut.chain_id+"_"+mut.mutation_event+".fasta"
    
    # downloading and cleaning PDB structure
    pdbdata.download_structure(pdb_id=mut.pdb_id)
    pdbdata.clean(pdb_id=mut.pdb_id, chain_id=mut.chain_id, selector=ChainAndAminoAcidSelect(mut.chain_id))

    # getting 0-based mutation site
    zero_based_mutation_site = Utils.get_zero_based_mutation_site(cln_pdb_file, mut.chain_id, mut.mutation_site)
    print("Row no:{}->{}{}, mutation:{}, 0-based_mutaiton_site:{}".format(i+1, mut.pdb_id, mut.chain_id, mut.mutation_event, zero_based_mutation_site))

    # generating wild and mutant fasta file
    pdbdata.generate_fasta_from_pdb(pdb_id=mut.pdb_id, chain_id=mut.chain_id, input_pdb_filepath=cln_pdb_file, save_as_fasta=True, output_fasta_dir=fastas_dir)
    pdbdata.create_mutant_fasta_file(wild_fasta_file=wild_fasta_file, mutant_fasta_file=mutant_fasta_file, zero_based_mutation_site=zero_based_mutation_site, mutant_residue=mut.mutant_residue, mutation=mut.mutation_event)
    
    print()
    del mut
    if i+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break

