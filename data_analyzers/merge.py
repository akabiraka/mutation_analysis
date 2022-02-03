import sys
sys.path.append("../mutation_analysis")
import pandas as pd
import os

class Mutation(object):
    def __init__(self) -> None:
        self.merged_file_path = "merged.csv"
        self.cols = ["pdb_id", "chain_id", "uniprot_id", "mutation_event", "wild_residue", "mutant_residue", "mutation_site", 
        "ddg", "ph", "temp", "method", "source_file_path", "source_id", "source_row_index", "pubmed_id", "extra_info"]
        self.merged_df = self.__init_merged_df()

        self.pdb_id = ""
        self.chain_id = ""
        self.mutation_event = ""
        self.wild_residue = ""
        self.mutant_residue = ""
        self.mutation_site = -1
        self.ddg = 0.0
        self.ph = 0.0
        self.temp = 0.0
        self.method = ""
        self.source_file_path = ""
        self.source_id = ""
        self.source_row_index = -1
        self.uniprot_id = ""
        self.pubmed_id = ""
        self.extra_info = ""

    def __str__(self):
        return "pdb_id={}, uniprot_id={}, chain_id={}, mutation_event={}, ddg={}, ph={}, temp(c)={}, method={}, source_file_path={}, source_id={}, source_row_index={}, pubmed_id={}, extra_info={}".format(
            self.pdb_id, self.uniprot_id, self.chain_id, self.mutation_event, self.ddg, self.ph, self.temp, 
            self.method, self.source_file_path, self.source_id, self.source_row_index, self.pubmed_id, self.extra_info)
    

    def __init_merged_df(self):
        if os.path.exists(self.merged_file_path): 
            self.merged_df = pd.read_csv(self.merged_file_path)
        else: 
            self.merged_df = pd.DataFrame(columns=self.cols)
            self.save()
        return self.merged_df


    def save(self):
        self.merged_df.to_csv(self.merged_file_path, index=False)


    def merge(self):
        # print(self)
        check_unique_mask = (self.merged_df["pdb_id"]==self.pdb_id) & (self.merged_df["chain_id"]==self.chain_id) & (self.merged_df["mutation_event"]==self.mutation_event) \
            & (self.merged_df["ddg"]==self.ddg) & (self.merged_df["ph"]==self.ph) & (self.merged_df["temp"]==self.temp) & (self.merged_df["method"]==self.method)
        if self.merged_df[check_unique_mask].size == 0:
            data_dict = {"pdb_id":self.pdb_id, "uniprot_id":self.uniprot_id, "chain_id":self.chain_id, "mutation_event":self.mutation_event, 
            "wild_residue": self.wild_residue, "mutant_residue": self.mutant_residue, "mutation_site": self.mutation_site, 
            "ddg": self.ddg, "ph": self.ph, "temp": self.temp, "method": self.method, "source_file_path": self.source_file_path, 
            "source_id": self.source_id, "source_row_index":self.source_row_index, "pubmed_id": self.pubmed_id, "extra_info": self.extra_info}
            self.merged_df = self.merged_df.append(data_dict, ignore_index=True)
            print(self.merged_df.shape)


class ProThermDB(object):
    def __init__(self) -> None:
        self.file_path = "data/ProThermDB_Oct_2021.tsv"
        self.df = pd.read_csv(self.file_path, delimiter="\t")

    def get_mutation(self, row):
        if row.ddg=="-": return
        mutation = Mutation()
        mutation.pdb_id = row.PDB_wild
        mutation.mutation_event = ""
        mutation.chain_id = ""
        mutation.ph = row.pH
        mutation.temp = row.T_C
        mutation.method = row.METHOD
        mutation.ddg = row.ddg
        mutation.extra_info = row.MUTATION+","+row.PDB_Chain_Mutation

        mutation.source_file_path = self.file_path
        mutation.source_id = row.NO
        mutation.source_row_index = row.Index
        mutation.uniprot_id = row.UniProt_ID
        mutation.pubmed_id = row.PubMed_ID
        return mutation



proThermDB = ProThermDB()
# print(proThermDB.df.head())
for row in proThermDB.df.itertuples():
    # print(row)
    mutation = proThermDB.get_mutation(row)
    if isinstance(mutation, Mutation): 
        mutation.merge()
        mutation.save()
    if row.Index == 185: break