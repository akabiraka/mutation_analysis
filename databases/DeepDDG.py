import sys
sys.path.append("../mutation_analysis")

import pandas as pd
import math
from objects.Mutation import Mutation
from databases.I_Database import I_Database

class DeepDDG(I_Database):
    def __init__(self, inp_file_path, out_file_path, sheet_name="sheet1") -> None:
        self.inp_file_path = inp_file_path
        df = pd.read_excel(self.inp_file_path, sheet_name=sheet_name)
        super().__init__(df, out_file_path)
        # print(self.df.head())

    def parse_mutation_event(self, mutation_event):
        mutation_event = mutation_event.split("_")
        return mutation_event[0], int(mutation_event[1]), mutation_event[2]

    def validate_pdb_id(self, pdb_id):
        return pdb_id[:4]

    def get_mutations(self, row):
        mutation = Mutation()
        mutation.pdb_id = self.validate_pdb_id(row.PDB_ID_with_modifications_to_be_made)
        mutation.chain_id = None
        mutation.uniprot_id = row.Uniprot_ID
        mutation.pubmed_id = row.PubMed_ID
        mutation.protein = row.Protein_name

        
        mutation.wild_residue, mutation.mutation_site, mutation.mutant_residue = self.parse_mutation_event(row.Mutation)
        mutation.mutation_event = mutation.wild_residue + str(mutation.mutation_site) + mutation.mutant_residue
    
        mutation.ddg = self.validate_ddg(row.ddg)
        mutation.dtm = None
        mutation.ph = self.validate_ph(row.pH)
        mutation.temp = self.validate_temp(row.T)

        mutation.inverse_pdb_id = None
        mutation.inverse_chain_id = None

        mutation.method = None
        mutation.event_based_on = None

        mutation.source_file_path = self.inp_file_path
        mutation.source_id = None
        mutation.source_row_index = row.Index
        
        mutation.extra_info = row.Source
        # print(mutation)
        return [mutation]


inp_file_path = "data/downloaded_as/DeepDDG_train_test.xlsx"
out_file_path = "data/clean_1/DeepDDG_S276.csv"
deepDDG_test = DeepDDG(inp_file_path, out_file_path, sheet_name="test")
deepDDG_test.run(0, 100000)