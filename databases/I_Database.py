import sys
from tkinter.tix import Tree
sys.path.append("../mutation_analysis")

from abc import ABC, abstractmethod
from objects.Mutation import Mutation

class I_Database(ABC):
    def __init__(self, df, out_file_path) -> None:
        self.df = df
        self.out_file_path = out_file_path


    def validate_chain(self, chain):
        return chain

    def convert_to_float(self, x):
        try:
            return float(x)
        except ValueError:
            if type(x) == str and "," in x:
                comma_idx = x.find(",")
                return float(x[:comma_idx]+"."+x[comma_idx+1:])

    def validate_ddg(self, ddg):
        return self.convert_to_float(ddg)

    def validate_dtm(self, dtm):
        return self.convert_to_float(dtm)

    def validate_ph(self, ph):
        return self.convert_to_float(ph)
    
    def validate_temp(self, temp):
        return self.convert_to_float(temp)

    def parse_mutation_event(self, mutation_event):
        wild_residue, mutation_site, mutant_residue = mutation_event[0], mutation_event[1:-1], mutation_event[-1]
        if mutation_site.isdigit():
            hetero_flag = " "
            insertion_code = " "
            mutation_site = int(mutation_site)    
        elif mutation_site[:-1].isdigit():
            hetero_flag = " "
            insertion_code = mutation_site[-1]
            mutation_site = int(mutation_site[:-1])
        elif mutation_site.startswith("H_") or mutation_site.startswith("W"):
            raise NotImplementedError()
        return wild_residue, hetero_flag, mutation_site, insertion_code, mutant_residue

    @DeprecationWarning
    def has_insertion_code(self, pdb, mutation_event):
        # based on PoPMuSiC_2, Saraboji
        if pdb=="1LVE" and mutation_event=="V27BL":
            return True
        elif pdb=="1LVE" and mutation_event=="L27CQ":
            return True
        elif pdb=="1LVE" and mutation_event=="L27CN":
            return True
        elif pdb == "1LVE" and mutation_event=="Y27DD":
            return True
        # based on ThermoMuteDB_single
        if pdb=="1C9O" and mutation_event=="L67AL":
            return True
        elif pdb=="1C9O" and mutation_event=="L66LA":
            return True
        elif pdb=="1IMQ" and mutation_event=="T27TG":
            return True
        # elif pdb=="'-" and mutation_event=="AVIG":
        #     return mutation_event[0]+"6"+mutation_event[-1]
        # elif pdb=="'-" and mutation_event=="AIIG":
        #     return mutation_event[0]+"2"+mutation_event[-1]
        else:
            return False

    def run(self, n_rows_to_skip, n_rows_to_evalutate):
        for row in self.df.itertuples():
            if row.Index+1 <= n_rows_to_skip: continue
            print("Row no: {}".format(row.Index))
            
            mutations = self.get_mutations(row)
            if mutations is not None:
                for mutation in mutations:
                    if isinstance(mutation, Mutation): 
                        mutation.save(self.out_file_path)
                        
            if row.Index+1 == n_rows_to_skip+n_rows_to_evalutate: 
                break


    def get_mutations(self, row):
        pass