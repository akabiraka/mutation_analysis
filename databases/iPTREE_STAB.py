import sys
sys.path.append("../mutation_analysis")

import pandas as pd
from objects.Mutation import Mutation
from databases.I_Database import I_Database

class iPTREE_STAB(I_Database):
    def __init__(self, inp_file_path, out_file_path) -> None:
        self.inp_file_path = inp_file_path
        df = pd.read_excel(self.inp_file_path)
        super().__init__(df, out_file_path)
        # print(self.df.head())

    def get_mutations(self, row):
        mutation = Mutation()
        mutation.pdb_id = row.PDB
        mutation.chain_id = None
        mutation.uniprot_id = None
        mutation.pubmed_id = None
        mutation.protein = None

        
        mutation.mutation_event = self.validate_mutation(mutation.pdb_id, row.Variation)
        mutation.wild_residue, mutation.mutation_site, mutation.mutant_residue = self.parse_mutation_event(mutation.mutation_event)
        
        mutation.ddg = self.validate_ddg(row.ddG)
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
        
        mutation.extra_info = None
        return [mutation]


inp_file_path = "data/downloaded_as/iPTREE-STAB_S1859.xlsx"
out_file_path = "data/clean_1/iPTREE_STAB.csv"

iptree_stab = iPTREE_STAB(inp_file_path, out_file_path)
iptree_stab.run(0, 100000)