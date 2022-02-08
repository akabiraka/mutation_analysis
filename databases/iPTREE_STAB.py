import sys
sys.path.append("../mutation_analysis")

import pandas as pd
from databases.Mutation import Mutation

class iPTREE_STAB(object):
    def __init__(self, file_path) -> None:
        self.file_path = file_path
        self.df = pd.read_excel(self.file_path)
        # print(self.df.head())

    def get_mutations(self, row):
        mutation = Mutation()
        mutation.pdb_id = row.PDB
        mutation.chain_id = None
        mutation.uniprot_id = None
        mutation.pubmed_id = None
        mutation.protein = None

        mutation_event = row.Variation
        mutation.mutation_event = mutation_event
        mutation.wild_residue = mutation_event[0]
        mutation.mutation_site = int(mutation_event[1:-1])
        mutation.mutant_residue = mutation_event[-1]
        
        mutation.ddg = float(row.ddG)
        mutation.ph = float(row.pH)
        mutation.temp = float(row.T)

        mutation.inverse_pdb_id = None
        mutation.inverse_chain_id = None

        mutation.method = None
        mutation.event_based_on = None

        mutation.source_file_path = self.file_path
        mutation.source_id = None
        mutation.source_row_index = row.Index
        
        mutation.extra_info = None
        return [mutation]


inp_file_path = "data/downloaded_as/iPTREE-STAB_S1859.xlsx"
out_file_path = "data/clean_1/iPTREE_STAB.csv"

iptree_stab = iPTREE_STAB(inp_file_path)
n_rows_to_skip = 0
n_rows_to_evalutate = 1000000

for row in iptree_stab.df.itertuples():
    if row.Index+1 <= n_rows_to_skip: continue
    print(row.Index)
    
    mutations = iptree_stab.get_mutations(row)
    if mutations is not None:
        for mutation in mutations:
            if isinstance(mutation, Mutation): 
                mutation.save(out_file_path)
                
    if row.Index+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break