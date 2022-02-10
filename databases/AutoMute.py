import sys
sys.path.append("../mutation_analysis")

import pandas as pd
from databases.Mutation import Mutation
from databases.I_Database import I_Database

class AutoMute(I_Database):
    def __init__(self, file_path) -> None:
        super().__init__()
        self.file_path = file_path
        self.df = pd.read_excel(self.file_path)
        # print(self.df.head())

    def get_mutations(self, row):
        mutation = Mutation()
        mutation.pdb_id = row.PDB
        mutation.chain_id = self.validate_chain(row.chain)
        mutation.uniprot_id = None
        mutation.pubmed_id = None
        mutation.protein = None

        mutation.mutation_event = self.validate_mutation(mutation.pdb_id, row.Variation)
        mutation.wild_residue, mutation.mutation_site, mutation.mutant_residue = self.parse_mutation_event(mutation.mutation_event)
    
        mutation.ddg = self.validate_ddg(row.ddG)
        mutation.dtm = None
        mutation.ph = self.validate_ph(row.pH)
        mutation.temp = self.validate_temp(row.temp)

        mutation.inverse_pdb_id = None
        mutation.inverse_chain_id = None

        mutation.method = None
        mutation.event_based_on = None

        mutation.source_file_path = self.file_path
        mutation.source_id = None
        mutation.source_row_index = row.Index
        
        mutation.extra_info = None
        # print(mutation)
        return [mutation]


inp_file_path = "data/downloaded_as/AUTOMUTE_S1925.xlsx"
out_file_path = "data/clean_1/AUTOMUTE_S1925.csv"

# inp_file_path = "data/downloaded_as/AUTOMUTE_S1962.xlsx"
# out_file_path = "data/clean_1/AUTOMUTE_S1962.csv"

automute = AutoMute(inp_file_path)
n_rows_to_skip = 0
n_rows_to_evalutate = 100#@0000

for row in automute.df.itertuples():
    if row.Index+1 <= n_rows_to_skip: continue
    print(row.Index)
    
    mutations = automute.get_mutations(row)
    if mutations is not None:
        for mutation in mutations:
            if isinstance(mutation, Mutation): 
                mutation.save(out_file_path)
                
    if row.Index+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break