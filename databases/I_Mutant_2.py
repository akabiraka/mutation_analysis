import sys
sys.path.append("../mutation_analysis")

import pandas as pd
from databases.Mutation import Mutation

class I_Mutant_2(object):
    def __init__(self, file_path) -> None:
        self.file_path = file_path
        self.df = pd.read_excel(self.file_path)
        # print(self.df.head())

    def get_mutations(self, row):
        mutation = Mutation()
        mutation.pdb_id = row.PDB
        mutation.mutation_event = row.Variation
        mutation.wild_residue = row.Variation[0]
        mutation.mutant_residue = row.Variation[-1]
        mutation.mutation_site = int(row.Variation[1:-1])
        mutation.ddg = float(row.ddG)
        mutation.ph = float(row.pH)
        mutation.temp = float(row.Temperature)
        mutation.source_file_path = self.file_path
        mutation.source_row_index = row.Index
        return [mutation]


i_Mutant_2_seq = I_Mutant_2("data/downloaded_as/I_Mutant2.0_seq_S2087.xlsx")
n_rows_to_skip = 0
n_rows_to_evalutate = 1000000
out_file_path="data/clean_1/I_Mutant_2_seq.csv"

for row in i_Mutant_2_seq.df.itertuples():
    if row.Index+1 <= n_rows_to_skip: continue
    print(row.Index)
    
    mutations = i_Mutant_2_seq.get_mutations(row)
    if mutations is not None:
        for mutation in mutations:
            if isinstance(mutation, Mutation): 
                mutation.save(out_file_path)
                
    
    if row.Index+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break

i_Mutant_2_struct = I_Mutant_2("data/downloaded_as/I_Mutant2.0_S1948_structure.xlsx")
n_rows_to_skip = 0
n_rows_to_evalutate = 1000000
out_file_path="data/clean_1/I_Mutant_2_structure.csv"

for row in i_Mutant_2_struct.df.itertuples():
    if row.Index+1 <= n_rows_to_skip: continue
    print(row.Index)
    
    mutations = i_Mutant_2_struct.get_mutations(row)
    if mutations is not None:
        for mutation in mutations:
            if isinstance(mutation, Mutation): 
                mutation.save(out_file_path)
                
    
    if row.Index+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break