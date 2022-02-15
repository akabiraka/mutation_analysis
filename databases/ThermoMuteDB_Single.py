import sys
sys.path.append("../mutation_analysis")

import pandas as pd
from databases.Mutation import Mutation
from databases.I_Database import I_Database

class ThermoMuteDB_Single(I_Database):
    def __init__(self, inp_file_path, out_file_path) -> None:
        self.inp_file_path = inp_file_path
        df = pd.read_csv(self.inp_file_path)
        super().__init__(df, out_file_path)
        # print(self.df.head())


    def validate_temp(self, temp):
        if temp=="'-": return
        else: return float(temp) - 273.15#converting to Celsius

    def validate_ddg(self, ddg):
        if type(ddg)==str and ddg[0]=="'": return self.convert_to_float(ddg[1:])
        else: return self.convert_to_float(ddg)

    def validate_dtm(self, dtm):
        if dtm=="'-": return float("nan")
        if type(dtm)==str and dtm[0]=="'": return self.convert_to_float(dtm[1:])
        else: return self.convert_to_float(dtm)


    def get_mutations(self, row):
        if row.ddG=="'-": return
        if row.PDB=="'-": return

        mutation = Mutation()
        mutation.pdb_id = row.PDB
        mutation.chain_id = None
        mutation.uniprot_id = None
        mutation.pubmed_id = row.Reference
        mutation.protein = row.Protein_Name

        mutation.mutation_event = self.validate_mutation(mutation.pdb_id, row.Mutation)
        mutation.wild_residue, mutation.mutation_site, mutation.mutant_residue = self.parse_mutation_event(mutation.mutation_event)
        
        mutation.ddg = self.validate_ddg(row.ddG)
        mutation.dtm = self.validate_dtm(row.dtm)
        mutation.ph = self.validate_ph(row.pH)
        mutation.temp = self.validate_temp(row.T)

        mutation.inverse_pdb_id = None
        mutation.inverse_chain_id = None

        mutation.method = row.Method if row.Method!="Unavailable" else None
        mutation.event_based_on = None

        mutation.source_file_path = self.inp_file_path
        mutation.source_id = row.ID
        mutation.source_row_index = row.Index
        
        mutation.extra_info = None
        return [mutation]


inp_file_path = "data/downloaded_as/ThermoMutDB_single_point.csv"
out_file_path = "data/clean_1/ThermoMutDB_single.csv"

thermoMuteDB_Single = ThermoMuteDB_Single(inp_file_path, out_file_path)
thermoMuteDB_Single.run(0, 100000)