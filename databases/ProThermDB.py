import sys
sys.path.append("../mutation_analysis")
import pandas as pd

from databases.Mutation import Mutation
from databases.I_Database import I_Database

class ProThermDB(I_Database):
    def __init__(self, inp_file_path, out_file_path) -> None:
        self.inp_file_path = inp_file_path
        df = pd.read_csv(self.inp_file_path, delimiter="\t")
        super().__init__(df, out_file_path)
        # print(self.df.head())
        
    def validate_ddg(self, ddg):
        if ddg == "-": return None
        return float(ddg.split(" ")[0]) #case: -0.76 (0.9)

    def validate_dtm(self, dtm):
        if dtm == "-": return None
        return float(dtm.split(" ")[0]) #case: 61.9 (0.4)

    def validate_ph(self, ph):
        if ph == "-": return None
        return self.convert_to_float(ph)

    def validate_temp(self, temp):
        if temp == "-": return None
        return self.convert_to_float(temp)
    

    def __populate_general_info(self, mutation, mutation_event, row):
        mutation.pdb_id = row.PDB_wild
        mutation.uniprot_id = row.UniProt_ID
        mutation.pubmed_id = row.PubMed_ID
        mutation.protein = row.PROTEIN

        mutation.wild_residue, mutation.mutation_site, mutation.mutant_residue = self.parse_mutation_event(mutation_event)
        mutation.mutation_event = mutation.wild_residue+str(mutation.mutation_site)+mutation.mutant_residue

        mutation.ddg = self.validate_ddg(row.ddg)
        mutation.dtm = self.validate_dtm(row.dtm)
        mutation.ph = self.validate_ph(row.pH)
        mutation.temp = self.validate_temp(row.T_C)

        mutation.inverse_pdb_id = None
        mutation.inverse_chain_id = None

        mutation.method = row.METHOD
        mutation.event_based_on = row.MUTATION.split("(")[1].rstrip()[8:-1]
        
        mutation.source_file_path = self.inp_file_path
        mutation.source_id = row.NO
        mutation.source_row_index = row.Index

        mutation.extra_info = row.MUTATION+","+row.PDB_Chain_Mutation
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
        
inp_file_path = "data/downloaded_as/ProThermDB_Oct_2021.tsv"
out_file_path="data/clean_1/ProThermDB.csv"
proThermDB = ProThermDB(inp_file_path, out_file_path)
proThermDB.run(0, 100000)