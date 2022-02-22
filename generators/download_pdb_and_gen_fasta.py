import sys

from numpy import mat
sys.path.append("../mutation_analysis")

import math
import pandas as pd
from objects.PDBData import PDBData
from objects.Selector import ChainAndAminoAcidSelect
import generators.utils as Utils
from objects.Mutation import Mutation

# configurations
pdb_dir = "data/pdbs/"
pdbs_clean_dir = "data/pdbs_clean/"
fastas_dir = "data/fastas/"
CIF = "mmCif"

inp_file_path = "data/clean_1/DeepDDG_S276.csv"
out_file_path = "data/clean_2/DeepDDG_S276.csv"
n_rows_to_skip = 0
n_rows_to_evalutate = 100000

# object initialization
pdbdata = PDBData(pdb_dir=pdb_dir)

df = pd.read_csv(inp_file_path)

def integrat_chain_id(pdb_id, chain_id):
    if type(chain_id)==str and len(chain_id)==1: return chain_id

    # these are in DeepDDG_test_276 set
    if pdb_id=="1gua": return "B"
    if pdb_id=="1iv7": return "B"
    if pdb_id=="2clr": return "B"

    chain_id = pdbdata.get_A_or_first_chain_id(pdb_id=pdb_id)
    return chain_id

for i, row in df.iterrows():
    if i+1 <= n_rows_to_skip: continue

    # getting mutation object
    mut = Mutation().get_mutation(row)

    # downloading PDB
    pdbdata.download_structure(pdb_id=mut.pdb_id)

    # updating chain_id
    mut.chain_id = integrat_chain_id(mut.pdb_id, mut.chain_id)
    df.loc[i, "chain_id"] = mut.chain_id
  
    # creating necessary file paths
    cln_pdb_file = pdbs_clean_dir+mut.pdb_id+mut.chain_id+".pdb"
    wild_fasta_file = fastas_dir+mut.pdb_id+mut.chain_id+".fasta"
    mutant_fasta_file = fastas_dir+mut.pdb_id+mut.chain_id+"_"+mut.mutation_event+".fasta"
    
    # cleaning PDB structure
    pdbdata.clean(pdb_id=mut.pdb_id, chain_id=mut.chain_id, selector=ChainAndAminoAcidSelect(mut.chain_id))

    # checking if mutation site has expected residue in PDB and specified chain
    flag, pdb_residue_name = pdbdata.does_mutation_site_has_expected_residue(cln_pdb_file, mut.chain_id, mut.mutation_site, mut.wild_residue)
    if not flag: 
        print(mut)
        print(pdb_residue_name)
        break
    
    zero_based_mutation_site = pdbdata.get_zero_based_mutation_site(cln_pdb_file, mut.chain_id, mut.mutation_site)
    print("Row no:{}->{}{}, mutation:{}, 0-based_mutaiton_site:{}".format(i+1, mut.pdb_id, mut.chain_id, mut.mutation_event, zero_based_mutation_site))

    # generating wild and mutant fasta file
    pdbdata.generate_fasta_from_pdb(pdb_id=mut.pdb_id, chain_id=mut.chain_id, input_pdb_filepath=cln_pdb_file, save_as_fasta=True, output_fasta_dir=fastas_dir)
    pdbdata.create_mutant_fasta_file(wild_fasta_file=wild_fasta_file, mutant_fasta_file=mutant_fasta_file, zero_based_mutation_site=zero_based_mutation_site, mutant_residue=mut.mutant_residue, mutation=mut.mutation_event)
    
    print()
    del mut
    if i+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break

df.to_csv(out_file_path, index=False)