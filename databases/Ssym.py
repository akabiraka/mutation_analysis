import sys
sys.path.append("../mutation_analysis")

import pandas as pd
from databases.Mutation import Mutation
from Bio.PDB.Polypeptide import three_to_one

class Ssym(object):
    def __init__(self, file_path) -> None:
        self.file_path = file_path
        self.df = pd.read_excel(self.file_path)
        # print(self.df.head())

    def __convert_to_float(self, x):
        try:
            return float(x)
        except ValueError:
            if type(x) == str and "," in x:
                comma_idx = x.find(",")
                return float(x[:comma_idx]+"."+x[comma_idx+1:])
        

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
        
        mutation.ddg = self.__convert_to_float(row.ddg)
        mutation.ph = self.__convert_to_float(row.pH)
        mutation.temp = self.__convert_to_float(row.T)

        mutation.inverse_pdb_id = row.inv_PDB
        mutation.inverse_chain_id = row.inv_Chain

        mutation.method = None
        mutation.event_based_on = None

        mutation.source_file_path = self.file_path
        mutation.source_id = None
        mutation.source_row_index = row.Index
        
        mutation.extra_info = None
        return [mutation]


ssym = Ssym("data/downloaded_as/Ssym.xlsx")
n_rows_to_skip = 0
n_rows_to_evalutate = 1000000
out_file_path="data/clean_1/Ssym.csv"

for row in ssym.df.itertuples():
    if row.Index+1 <= n_rows_to_skip: continue
    print(row.Index)
    
    mutations = ssym.get_mutations(row)
    if mutations is not None:
        for mutation in mutations:
            if isinstance(mutation, Mutation): 
                mutation.save(out_file_path)
                
    
    if row.Index+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break