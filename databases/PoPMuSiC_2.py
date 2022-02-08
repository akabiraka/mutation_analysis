import sys
sys.path.append("../mutation_analysis")

import pandas as pd
from databases.Mutation import Mutation

class PoPMuSiC_2(object):
    def __init__(self, file_path) -> None:
        self.file_path = file_path
        self.df = pd.read_excel(self.file_path)
        # print(self.df.head())

    def get_mutations(self, row):
        if row.Index in range(1315, 1319): return #skipped
        mutation = Mutation()
        mutation.pdb_id = row.PDB[:-1]
        mutation.chain_id = row.CHAIN
        mutation.mutation_event = row.Variation
        mutation.event_based_on = None
        mutation.wild_residue = row.Variation[0]
        mutation.mutant_residue = row.Variation[-1]
        mutation.mutation_site = int(row.Variation[1:-1])
        mutation.ddg = float(row.ddG)
        mutation.ph = float(row.pH)
        mutation.temp = float(row.T)
        mutation.method = None
        mutation.source_file_path = self.file_path
        mutation.source_id = None
        mutation.source_row_index = row.Index
        mutation.uniprot_id = None
        mutation.pubmed_id = None
        mutation.extra_info = None
        mutation.protein = None
        return [mutation]


poPMuSiC_2 = PoPMuSiC_2("data/downloaded_as/PoPMuSiC-2.0_S2648.xlsx")
n_rows_to_skip = 0
n_rows_to_evalutate = 1000000
out_file_path="data/clean_1/PoPMuSiC_2.csv"

for row in poPMuSiC_2.df.itertuples():
    if row.Index+1 <= n_rows_to_skip: continue
    print(row.Index)
    
    mutations = poPMuSiC_2.get_mutations(row)
    if mutations is not None:
        for mutation in mutations:
            if isinstance(mutation, Mutation): 
                mutation.save(out_file_path)
                
    
    if row.Index+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break