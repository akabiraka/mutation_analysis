import sys
sys.path.append("../mutation_analysis")

import pandas as pd
from databases.Mutation import Mutation

class FireProtDB(object):
    def __init__(self, file_path) -> None:
        self.file_path = file_path
        self.df = pd.read_csv(self.file_path)
        # print(self.df.head())

    def get_mutations(self, row):
        mutation = Mutation()
        mutation.pdb_id = row.pdb_id
        mutation.chain_id = row.chain
        mutation.uniprot_id = row.uniprot_id
        mutation.pubmed_id = row.publication_pubmed
        mutation.protein = row.protein_name

        mutation.wild_residue = row.wild_type
        mutation.mutation_site = int(row.position)
        mutation.mutant_residue = row.mutation
        mutation.mutation_event = mutation.wild_residue + str(mutation.mutation_site) + mutation.mutant_residue
        
        mutation.ddg = float(row.ddG)
        mutation.dtm = float(row.dTm)
        mutation.ph = float(row.pH)
        mutation.temp = float(row.tm)

        mutation.inverse_pdb_id = None
        mutation.inverse_chain_id = None

        mutation.method = row.method
        mutation.event_based_on = None

        mutation.source_file_path = self.file_path
        mutation.source_id = row.experiment_id
        mutation.source_row_index = row.Index
        
        mutation.extra_info = None
        return [mutation]


inp_file_path = "data/downloaded_as/FireProtDB.csv"
out_file_path = "data/clean_1/FireProtDB.csv"

fireProtDB = FireProtDB(inp_file_path)
n_rows_to_skip = 0
n_rows_to_evalutate = 1000000

for row in fireProtDB.df.itertuples():
    if row.Index+1 <= n_rows_to_skip: continue
    print(row.Index)
    
    mutations = fireProtDB.get_mutations(row)
    if mutations is not None:
        for mutation in mutations:
            if isinstance(mutation, Mutation): 
                mutation.save(out_file_path)
                
    if row.Index+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break