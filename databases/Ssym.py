import sys
sys.path.append("../mutation_analysis")

import pandas as pd
from objects.Mutation import Mutation
from Bio.PDB.Polypeptide import three_to_one
from databases.I_Database import I_Database

class Ssym(I_Database):
    def __init__(self, inp_file_path, out_file_path) -> None:
        self.inp_file_path = inp_file_path
        df = pd.read_excel(self.inp_file_path)
        super().__init__(df, out_file_path)
        

    def get_mutations(self, row):
        mutation = Mutation()
        mutation.pdb_id = row.PDB
        mutation.chain_id = row.Chain
        mutation.uniprot_id = None
        mutation.pubmed_id = row.PubMed_ID
        mutation.protein = None

        
        mutation.wild_residue = three_to_one(row.Wild_Type)
        mutation.mutation_site = int(row.Res_Num)
        mutation.mutant_residue = three_to_one(row.Mutant)
        mutation.mutation_event = mutation.wild_residue+str(mutation.mutation_site)+mutation.mutant_residue
        
        mutation.ddg = self.validate_ddg(row.ddg)
        mutation.dtm = None
        mutation.ph = self.validate_ph(row.pH)
        mutation.temp = self.validate_temp(row.T)

        mutation.inverse_pdb_id = row.inv_PDB
        mutation.inverse_chain_id = row.inv_Chain

        mutation.method = None
        mutation.event_based_on = None

        mutation.source_file_path = self.inp_file_path
        mutation.source_id = None
        mutation.source_row_index = row.Index
        
        mutation.extra_info = None
        return [mutation]

inp_file_path = "data/downloaded_as/Ssym.xlsx"
out_file_path="data/clean_1/Ssym.csv"
ssym = Ssym(inp_file_path, out_file_path)
ssym.run(0, 10000)