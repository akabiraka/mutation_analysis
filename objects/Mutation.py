import sys
sys.path.append("../mutation_analysis")
import pandas as pd
import os
import math

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

    def get_mutation(self, row):
        self.pdb_id = row.pdb_id.lower() if type(row.pdb_id)==str else row.pdb_id
        self.chain_id = row.chain_id.upper() if type(row.chain_id)==str else row.chain_id
        self.uniprot_id = row.uniprot_id
        self.pubmed_id = row.pubmed_id
        self.protein = row.protein

        self.mutation_event = row.mutation_event
        self.wild_residue = row.wild_residue
        self.mutation_site = row.mutation_site
        self.mutant_residue = row.mutant_residue

        self.ddg = row.ddg
        self.dtm = row.dtm
        self.ph = row.ph
        self.temp = row.temp

        self.inverse_pdb_id = row.inverse_pdb_id.upper() if type(row.inverse_pdb_id)==str else row.inverse_pdb_id
        self.inverse_chain_id = row.inverse_chain_id.upper() if type(row.inverse_chain_id)==str else row.inverse_chain_id

        self.method = row.method
        self.event_based_on = row.event_based_on

        self.source_file_path = row.source_file_path
        self.source_id = row.source_id
        self.source_row_index = row.source_row_index
        
        self.extra_info = row.extra_info
        return self

    def __get_data_dict(self):
        return {"pdb_id":self.pdb_id.lower(), 
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

    def __should_update(self, new_value, prev_value):
        x = new_value!=prev_value and type(new_value)==str and type(prev_value)!=str and math.isnan(prev_value)
        if x: print(new_value, prev_value)
        return x

    def __update_current_row(self, out_df, idx, out_file_path):
        update_flag = False
        if self.__should_update(self.chain_id, out_df.loc[idx, "chain_id"]):
            out_df.loc[idx, "chain_id"] = self.chain_id
            update_flag = True
            
        if self.__should_update(self.uniprot_id, out_df.loc[idx, "uniprot_id"]):
            out_df.loc[idx, "uniprot_id"] = self.uniprot_id
            update_flag = True
            
        if self.__should_update(self.pubmed_id, out_df.loc[idx, "pubmed_id"]):
            out_df.loc[idx, "pubmed_id"] = self.pubmed_id
            update_flag = True
        
        if self.__should_update(self.protein, out_df.loc[idx, "protein"]):
            out_df.loc[idx, "protein"] = self.protein
            update_flag = True
        
        dtm = out_df.loc[idx, "dtm"]
        if self.dtm!=dtm and self.dtm!=None and not math.isnan(self.dtm) and math.isnan(dtm):
            out_df.loc[idx, "dtm"] = self.dtm
            update_flag = True
            
        if self.__should_update(self.inverse_pdb_id, out_df.loc[idx, "inverse_pdb_id"]):
            out_df.loc[idx, "inverse_pdb_id"] = self.inverse_pdb_id
            update_flag = True
        
        if self.__should_update(self.inverse_chain_id, out_df.loc[idx, "inverse_chain_id"]):
            out_df.loc[idx, "inverse_chain_id"] = self.inverse_chain_id
            update_flag = True
        
        if self.__should_update(self.method, out_df.loc[idx, "method"]):
            out_df.loc[idx, "method"] = self.method
            update_flag = True
        
        if self.__should_update(self.event_based_on, out_df.loc[idx, "event_based_on"]):
            out_df.loc[idx, "event_based_on"] = self.event_based_on
            update_flag = True
            
        
        if update_flag:
            out_df.loc[idx, "source_file_path"] = out_df.loc[idx, "source_file_path"] + ", " + self.source_file_path
            out_df.loc[idx, "source_id"] = str(out_df.loc[idx, "source_id"]) + ", " + str(self.source_id)
            out_df.loc[idx, "source_row_index"] = str(out_df.loc[idx, "source_row_index"]) + ", " + str(self.source_row_index)
            out_df.to_csv(out_file_path, index=False)
        return update_flag

    def __append_new_row(self, out_df, out_file_path):
        out_df = out_df.append(self.__get_data_dict(), ignore_index=True)
        out_df.to_csv(out_file_path, index=False)
        return True

    def save(self, out_file_path):
        out_df = self.__get_out_df(out_file_path)
        if (self.temp==None or math.isnan(self.temp)) and (self.ph==None or math.isnan(self.ph)):
            check_unique_mask = (out_df["pdb_id"]==self.pdb_id) & (out_df["mutation_event"]==self.mutation_event) & (out_df["ddg"]==self.ddg)
        elif self.temp==None or math.isnan(self.temp):
            check_unique_mask = (out_df["pdb_id"]==self.pdb_id) & (out_df["mutation_event"]==self.mutation_event) & (out_df["ddg"]==self.ddg) & (out_df["ph"]==self.ph)
        elif self.ph==None or math.isnan(self.ph): 
            check_unique_mask = (out_df["pdb_id"]==self.pdb_id) & (out_df["mutation_event"]==self.mutation_event) & (out_df["ddg"]==self.ddg) & (out_df["temp"]==self.temp)
        else: 
            check_unique_mask = (out_df["pdb_id"]==self.pdb_id) & (out_df["mutation_event"]==self.mutation_event) & (out_df["ddg"]==self.ddg) & (out_df["ph"]==self.ph) & (out_df["temp"]==self.temp)
        
        q_df = out_df[check_unique_mask]
        is_updated, is_appended = False, False
        if q_df.shape[0] == 1: #update existing row or not
            is_updated = self.__update_current_row(out_df, q_df.index[0], out_file_path)
        elif q_df.shape[0] == 0: #append new row
            is_appended = self.__append_new_row(out_df, out_file_path)

        return is_updated, is_appended

