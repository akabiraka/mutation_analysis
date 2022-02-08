import sys
sys.path.append("../mutation_analysis")

import pandas as pd
from databases.Mutation import Mutation

class PoPMuSiC_2(object):
    def __init__(self, file_path) -> None:
        self.file_path = file_path
        self.df = pd.read_excel(self.file_path)
        # print(self.df.head())

    def __validate_mutation(self, pdb, mutation_event):
        if pdb=="1LVE" and "27B" in mutation_event:
            return mutation_event[0]+"29"+mutation_event[-1]
        elif pdb=="1LVE" and "27C" in mutation_event:
            return mutation_event[0]+"30"+mutation_event[-1]
        elif pdb == "1LVE" and "27D" in mutation_event:
            return mutation_event[0]+"31"+mutation_event[-1]
        else:
            return mutation_event


    def get_mutations(self, row):
        mutation = Mutation()
        mutation.pdb_id = row.PDB[:-1]
        mutation.chain_id = row.CHAIN

        mutation_event = self.__validate_mutation(mutation.pdb_id, row.Variation)
        mutation.mutation_event = mutation_event
        mutation.wild_residue = mutation_event[0]
        mutation.mutation_site = int(mutation_event[1:-1])
        mutation.mutant_residue = mutation_event[-1]

        mutation.ddg = float(row.ddG)
        mutation.ph = float(row.pH)
        mutation.temp = float(row.T)
        mutation.source_file_path = self.file_path
        mutation.source_row_index = row.Index
        return [mutation]


poPMuSiC_2 = PoPMuSiC_2("data/downloaded_as/PoPMuSiC-2.0_S2648.xlsx")
n_rows_to_skip = 0
n_rows_to_evalutate = 1000000
out_file_path="data/clean_1/PoPMuSiC_2.csv"

for row in poPMuSiC_2.df.itertuples():
    if row.Index+1 <= n_rows_to_skip: continue
    print(row.Index)
    
    mutations = poPMuSiC_2.get_mutations(row)
    if mutations is not None:
        for mutation in mutations:
            if isinstance(mutation, Mutation): 
                mutation.save(out_file_path)
                
    
    if row.Index+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break