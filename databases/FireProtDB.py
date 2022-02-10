import sys
sys.path.append("../mutation_analysis")

import pandas as pd
from databases.Mutation import Mutation
from databases.I_Database import I_Database

class FireProtDB(I_Database):
    def __init__(self, inp_file_path, out_file_path) -> None:
        self.inp_file_path = inp_file_path
        df = pd.read_csv(self.inp_file_path)
        super().__init__(df, out_file_path)
        # print(self.df.head())


    def get_mutations(self, row):
        mutation = Mutation()
        mutation.pdb_id = row.pdb_id
        mutation.chain_id = self.validate_chain(row.chain)
        mutation.uniprot_id = row.uniprot_id
        mutation.pubmed_id = row.publication_pubmed
        mutation.protein = row.protein_name

        mutation.wild_residue = row.wild_type
        mutation.mutation_site = int(row.position)
        mutation.mutant_residue = row.mutation
        mutation.mutation_event = mutation.wild_residue + str(mutation.mutation_site) + mutation.mutant_residue
        
        mutation.ddg = self.validate_ddg(row.ddG)
        mutation.dtm = self.validate_dtm(row.dTm)
        mutation.ph = self.validate_ph(row.pH)
        mutation.temp = self.validate_temp(row.tm)
        
        mutation.inverse_pdb_id = None
        mutation.inverse_chain_id = None

        mutation.method = row.method
        mutation.event_based_on = None

        mutation.source_file_path = self.inp_file_path
        mutation.source_id = row.experiment_id
        mutation.source_row_index = row.Index
        
        mutation.extra_info = None
        # print(mutation)
        return [mutation]


inp_file_path = "data/downloaded_as/FireProtDB.csv"
out_file_path = "data/clean_1/FireProtDB.csv"

fireProtDB = FireProtDB(inp_file_path, out_file_path)
fireProtDB.run(0, 1000)
