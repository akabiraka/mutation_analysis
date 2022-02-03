import sys
sys.path.append("../mutation_analysis")
import pandas as pd
import os

class Mutation(object):
    def __init__(self) -> None:
        self.cols = ["pdb_id", "chain_id", "uniprot_id", "mutation_event", "event_based_on", "wild_residue", "mutant_residue", "mutation_site", 
                    "ddg", "ph", "temp", "method", "source_file_path", "source_id", "source_row_index", "pubmed_id", "extra_info"]
        
        self.pdb_id = None
        self.chain_id = None
        self.mutation_event = None
        self.event_based_on = None
        self.wild_residue = None
        self.mutant_residue = None
        self.mutation_site = None
        self.ddg = None
        self.ph = None
        self.temp = None
        self.method = None
        self.source_file_path = None
        self.source_id = None
        self.source_row_index = None
        self.uniprot_id = None
        self.pubmed_id = None
        self.extra_info = None


    def __str__(self):
        return "pdb_id={}, uniprot_id={}, chain_id={}, mutation_event={}, event_based_on={}, ddg={}, ph={}, temp(c)={}, method={}, source_file_path={}, source_id={}, source_row_index={}, pubmed_id={}, extra_info={}".format(
            self.pdb_id, self.uniprot_id, self.chain_id, self.mutation_event, self.event_based_on, self.ddg, self.ph, self.temp, 
            self.method, self.source_file_path, self.source_id, self.source_row_index, self.pubmed_id, self.extra_info)
    

    def __get_out_df(self, out_file_path):
        if os.path.exists(out_file_path): return pd.read_csv(out_file_path)
        return pd.DataFrame(columns=self.cols)
        

    def save(self, out_file_path):
        out_df = self.__get_out_df(out_file_path)
        check_unique_mask = (out_df["pdb_id"]==self.pdb_id) & (out_df["chain_id"]==self.chain_id) \
                            & (out_df["mutation_event"]==self.mutation_event) & (out_df["event_based_on"]==self.event_based_on) \
                            & (out_df["ddg"]==self.ddg) & (out_df["ph"]==self.ph) & (out_df["temp"]==self.temp) & (out_df["method"]==self.method)

        if out_df[check_unique_mask].size == 0:
            data_dict = {"pdb_id":self.pdb_id, 
                        "uniprot_id":self.uniprot_id, 
                        "chain_id":self.chain_id, 
                        "mutation_event":self.mutation_event, 
                        "event_based_on": self.event_based_on, 
                        "wild_residue": self.wild_residue, 
                        "mutant_residue": self.mutant_residue, 
                        "mutation_site": self.mutation_site, 
                        "ddg": self.ddg, 
                        "ph": self.ph, 
                        "temp": self.temp, 
                        "method": self.method, 
                        "source_file_path": self.source_file_path, 
                        "source_id": self.source_id, 
                        "source_row_index":self.source_row_index, 
                        "pubmed_id": self.pubmed_id, 
                        "extra_info": self.extra_info}
            
            out_df = out_df.append(data_dict, ignore_index=True)
            out_df.to_csv(out_file_path, index=False)
            print(out_df.shape)

