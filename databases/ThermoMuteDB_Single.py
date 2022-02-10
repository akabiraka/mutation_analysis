import sys
sys.path.append("../mutation_analysis")

import pandas as pd
from databases.Mutation import Mutation

class ThermoMuteDB_Single(object):
    def __init__(self, file_path) -> None:
        self.file_path = file_path
        self.df = pd.read_csv(self.file_path)
        # print(self.df.head())

    def __validate_ddg(self, ddg):
        if ddg=="'-": return
        if ddg[0]=="'": return float(ddg[1:])
        else: return float(ddg)

    def __validate_temp(self, temp):
        if temp=="'-": return
        else: return float(temp) - 273.15#converting to Celsius

    def __validate_mutation(self, pdb, mutation_event):
        if pdb=="1C9O" and mutation_event=="L67AL":
            return mutation_event[0]+"66"+mutation_event[-1]
        elif pdb=="1C9O" and mutation_event=="L66LA":
            return mutation_event[0]+"66"+mutation_event[-1]
        elif pdb=="1IMQ" and mutation_event=="T27TG":
            return mutation_event[0]+"27"+mutation_event[-1]
        elif pdb=="'-" and mutation_event=="AVIG":
            return mutation_event[0]+"6"+mutation_event[-1]
        elif pdb=="'-" and mutation_event=="AIIG":
            return mutation_event[0]+"2"+mutation_event[-1]
        else:
            return mutation_event


    def get_mutations(self, row):
        mutation = Mutation()
        mutation.pdb_id = row.PDB
        mutation.chain_id = None
        mutation.uniprot_id = None
        mutation.pubmed_id = row.Reference
        mutation.protein = row.Protein_Name

        mutation_event = self.__validate_mutation(mutation.pdb_id, row.Mutation)
        mutation.mutation_event = mutation_event
        mutation.wild_residue = mutation_event[0]
        mutation.mutation_site = int(mutation_event[1:-1])
        mutation.mutant_residue = mutation_event[-1]
        
        ddg = self.__validate_ddg(row.ddG)
        if ddg == None: return
        mutation.ddg = ddg
        mutation.ph = float(row.pH)
        mutation.temp = self.__validate_temp(row.T)

        mutation.inverse_pdb_id = None
        mutation.inverse_chain_id = None

        mutation.method = row.Method if row.Method!="Unavailable" else None
        mutation.event_based_on = None

        mutation.source_file_path = self.file_path
        mutation.source_id = row.ID
        mutation.source_row_index = row.Index
        
        mutation.extra_info = None
        return [mutation]


inp_file_path = "data/downloaded_as/ThermoMutDB_single_point.csv"
out_file_path = "data/clean_1/ThermoMutDB_single.csv"

thermoMuteDB_Single = ThermoMuteDB_Single(inp_file_path)
n_rows_to_skip = 0
n_rows_to_evalutate = 1000000

for row in thermoMuteDB_Single.df.itertuples():
    if row.Index+1 <= n_rows_to_skip: continue
    print(row.Index)
    
    mutations = thermoMuteDB_Single.get_mutations(row)
    if mutations is not None:
        for mutation in mutations:
            if isinstance(mutation, Mutation): 
                mutation.save(out_file_path)
                
    if row.Index+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break