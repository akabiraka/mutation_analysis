import sys
sys.path.append("../mutation_analysis")
import pandas as pd
import os

class Mutation(object):
    def __init__(self) -> None:
        self.cols = ["pdb_id", "chain_id", "uniprot_id", "pubmed_id", "protein",
                    "mutation_event", "wild_residue", "mutation_site", "mutant_residue",
                    "ddg", "dtm", "ph", "temp", 
                    "inverse_pdb_id", "inverse_chain_id",
                    "method", "event_based_on", 
                    "source_file_path", "source_id", "source_row_index", 
                    "extra_info"]
        
        self.pdb_id = None
        self.chain_id = None
        self.uniprot_id = None
        self.pubmed_id = None
        self.protein = None

        self.mutation_event = None
        self.wild_residue = None
        self.mutation_site = None
        self.mutant_residue = None
        
        self.ddg = None
        self.dtm = None
        self.ph = None
        self.temp = None

        self.inverse_pdb_id = None
        self.inverse_chain_id = None

        self.method = None
        self.event_based_on = None

        self.source_file_path = None
        self.source_id = None
        self.source_row_index = None
        
        self.extra_info = None

    def __get_data_dict(self):
        return {"pdb_id":self.pdb_id, 
                "chain_id":self.chain_id,
                "uniprot_id":self.uniprot_id, 
                "pubmed_id": self.pubmed_id, 
                "protein": self.protein,
                "mutation_event":self.mutation_event, 
                "wild_residue": self.wild_residue, 
                "mutation_site": self.mutation_site, 
                "mutant_residue": self.mutant_residue,
                "ddg": self.ddg, 
                "dtm": self.dtm,
                "ph": self.ph, 
                "temp": self.temp, 
                "inverse_pdb_id": self.inverse_pdb_id,
                "inverse_chain_id": self.inverse_chain_id,
                "method": self.method, 
                "event_based_on": self.event_based_on,
                "source_file_path": self.source_file_path, 
                "source_id": self.source_id, 
                "source_row_index":self.source_row_index, 
                "extra_info": self.extra_info}
                

    def __str__(self):
        return str(self.__get_data_dict())
    

    def __get_out_df(self, out_file_path):
        if os.path.exists(out_file_path): return pd.read_csv(out_file_path)
        return pd.DataFrame(columns=self.cols)
        

    def save(self, out_file_path):
        out_df = self.__get_out_df(out_file_path)
        check_unique_mask = (out_df["pdb_id"]==self.pdb_id) & (out_df["chain_id"]==self.chain_id) \
                            & (out_df["mutation_event"]==self.mutation_event) & (out_df["event_based_on"]==self.event_based_on) \
                            & (out_df["ddg"]==self.ddg) & (out_df["ph"]==self.ph) \
                            & (out_df["temp"]==self.temp) & (out_df["method"]==self.method)

        if out_df[check_unique_mask].size == 0:
            data_dict = self.__get_data_dict()
            # print(data_dict)
            out_df = out_df.append(data_dict, ignore_index=True)
            out_df.to_csv(out_file_path, index=False)
            # print(out_df.shape)

