import sys
sys.path.append("../mutation_analysis")

import pandas as pd
from databases.Mutation import Mutation

class PON_TStab():
    def __init__(self) -> None:
        self.file_path = "data/downloaded_as/PON-TStab_dataset.xlsx"
        self.df = pd.read_excel(self.file_path)


    def get_mutations(self, row):
        mutation = Mutation()
        mutation.pdb_id = row.PDB_in_Protherm
        mutation.chain_id = None
        mutation.mutation_event = row.Variation
        mutation.event_based_on = None
        mutation.wild_residue = row.Variation[0]
        mutation.mutant_residue = row.Variation[-1]
        mutation.mutation_site = int(row.Variation[1:-1])
        mutation.ddg = row.ddG
        mutation.ph = row.pH
        mutation.temp = row.T
        mutation.method = None
        mutation.source_file_path = self.file_path
        mutation.source_id = row.Record_id
        mutation.source_row_index = row.Index
        mutation.uniprot_id = None
        mutation.pubmed_id = row.Pubmed_id
        mutation.extra_info = None
        mutation.protein = row.protein_name
        return [mutation]


ponTStab = PON_TStab()
n_rows_to_skip = 0
n_rows_to_evalutate = 1000000
out_file_path="data/clean_1/PON_TStab.csv"

for row in ponTStab.df.itertuples():
    if row.Index+1 <= n_rows_to_skip: continue
    print(row.Index)
    
    mutations = ponTStab.get_mutations(row)
    if mutations is not None:
        for mutation in mutations:
            if isinstance(mutation, Mutation): 
                mutation.save(out_file_path)
                
    
    if row.Index+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break