import sys
sys.path.append("../mutation_analysis")
import math
import pandas as pd
from objects.Mutation import Mutation
from databases.I_Database import I_Database

class PoPMuSiC_2(I_Database):
    def __init__(self, inp_file_path, out_file_path) -> None:
        self.inp_file_path = inp_file_path
        df = pd.read_excel(self.inp_file_path)
        super().__init__(df, out_file_path)
        # print(self.df.head())


    def get_mutations(self, row):
        if math.isnan(row.ddG): return 

        mutation = Mutation()
        mutation.pdb_id = str(row.PDB[:-1])
        mutation.chain_id = self.validate_chain(row.CHAIN)
        mutation.uniprot_id = None
        mutation.pubmed_id = None
        mutation.protein = None

        
        mutation.mutation_event = row.Variation
        # mutation.wild_residue, mutation.mutation_site, mutation.mutant_residue = self.parse_mutation_event(mutation.mutation_event, self.has_insertion_code(mutation.pdb_id, row.Variation))
        mutation.wild_residue, mutation.hetero_flag, mutation.mutation_site, mutation.insertion_code, mutation.mutant_residue = self.parse_mutation_event(mutation.mutation_event)

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
        # print(mutation)
        return [mutation]

inp_file_path = "data/downloaded_as/PoPMuSiC-2.0_S2648.xlsx"
out_file_path = "data/clean_1/PoPMuSiC_2.csv"
poPMuSiC_2 = PoPMuSiC_2(inp_file_path, out_file_path)
poPMuSiC_2.run(0, 100000)