import sys
sys.path.append("../mutation_analysis")
import pandas as pd

from data_analyzers.Mutation import Mutation

class ProThermDB(object):
    def __init__(self) -> None:
        self.file_path = "data/downloaded_as/ProThermDB_Oct_2021.tsv"
        self.df = pd.read_csv(self.file_path, delimiter="\t")
        

    def __populate_general_info(self, mutation, mutation_event, row):
        mutation.wild_residue = mutation_event[0]
        mutation.mutant_residue = mutation_event[-1]
        mutation.mutation_site = int(mutation_event[1:-1])
        mutation.mutation_event = mutation.wild_residue+"_"+str(mutation.mutation_site)+"_"+mutation.mutant_residue

        mutation.event_based_on = row.MUTATION.split("(")[1].rstrip()[8:-1]
        mutation.method = row.METHOD
        mutation.pdb_id = row.PDB_wild
        if row.pH != "-": mutation.ph = float(row.pH)
        if row.T_C != "-": mutation.temp = float(row.T_C)
        
        mutation.ddg = float(row.ddg.split(" ")[0]) #case: -0.76 (0.9)
        
        mutation.extra_info = row.MUTATION+","+row.PDB_Chain_Mutation

        mutation.source_file_path = self.file_path
        mutation.source_id = row.NO
        mutation.source_row_index = row.Index
        mutation.uniprot_id = row.UniProt_ID
        mutation.pubmed_id = row.PubMed_ID

        return mutation


    def get_mutations(self, row):
        """Returns all mutations corresponding to a row.
        """
        if row.ddg=="-": return
        if row.PDB_wild=="-" and row.UniProt_ID=="-": return
        
        mutations = []
        if row.PDB_Chain_Mutation == "-":
            # cases:
            # row.MUTATION == D794N (Based on UniProt)
            # row.MUTATION == P94R,D160N,D165A,K190E,K196T(Based on Uniprot)
            # row.MUTATION == Q51K E84A (Based on UniProt)
            
            mutation_events = row.MUTATION.split("(")[0].rstrip()
            if "," in mutation_events: #P94R,D160N,D165A,K190E,K196T(Based on Uniprot)
                mutation_events = mutation_events.split(",")
            else: #Q51K E84A (Based on UniProt)
                mutation_events = mutation_events.split(" ")

            for mutation_event in mutation_events:
                mutation = Mutation()
                mutation = self.__populate_general_info(mutation, mutation_event, row)
                mutations.append(mutation)

        else:
            # cases:
            # row.PDB_Chain_Mutation == 1wq5_A:P28S
            # row.PDB_Chain_Mutation == 2lzm_A:W138Y 2lzm_A:W126Y 2lzm_A:W158Y
            # row.PDB_Chain_Mutation == 1yri:A_L42A
            # row.PDB_Chain_Mutation == 1ten_A:Q808K, 1ten_A:L820K, 1ten_A:D850K, 1ten_A:T890K
            
            if "," in row.PDB_Chain_Mutation: #1ten_A:Q808K, 1ten_A:L820K, 1ten_A:D850K, 1ten_A:T890K
                mutation_events = row.PDB_Chain_Mutation.split(",")
            else: #2lzm_A:W138Y 2lzm_A:W126Y 2lzm_A:W158Y
                mutation_events = row.PDB_Chain_Mutation.split(" ")

            for mutation_event in mutation_events:
                mutation_event = mutation_event.rstrip()
            
                mutation = Mutation()
                colon_idx, underscore_idx = mutation_event.find(":"), mutation_event.find("_")
                if colon_idx < underscore_idx: #1yri:A_L42A
                    mutation.chain_id = mutation_event[colon_idx+1:underscore_idx]
                    mutation_event = mutation_event[underscore_idx+1:]
                else: #1yri_A:L42A
                    mutation.chain_id = mutation_event[underscore_idx+1:colon_idx]
                    mutation_event = mutation_event[colon_idx+1:]

                mutation = self.__populate_general_info(mutation, mutation_event, row)
                mutations.append(mutation)
        
        return mutations
        

proThermDB = ProThermDB()
n_rows_to_skip = 0
n_rows_to_evalutate = 1000000
out_file_path="data/clean_1/ProThermDB.csv"

for row in proThermDB.df.itertuples():
    if row.Index+1 <= n_rows_to_skip: continue
    # print(row.Index)
    
    mutations = proThermDB.get_mutations(row)
    if mutations is not None:
        for mutation in mutations:
            if isinstance(mutation, Mutation): 
                mutation.save(out_file_path)
                
    
    if row.Index+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break