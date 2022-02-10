import sys
sys.path.append("../mutation_analysis")

from abc import ABC, abstractmethod

class I_Database(ABC):

    def init(self) -> None:
        super().init()

    def validate_chain(self, chain):
        # from AutoMute
        if chain=="@": return None
        else: return chain

    def convert_to_float(self, x):
        try:
            return float(x)
        except ValueError:
            if type(x) == str and "," in x:
                comma_idx = x.find(",")
                return float(x[:comma_idx]+"."+x[comma_idx+1:])

    def validate_ddg(self, ddg):
        if ddg=="'-": return float("nan")
        if type(ddg)==str and ddg[0]=="'": return self.convert_to_float(ddg[1:])
        else: return self.convert_to_float(ddg)

    def validate_dtm(self, dtm):
        return self.convert_to_float(dtm)


    def validate_ph(self, ph):
        return self.convert_to_float(ph)

    
    def validate_temp(self, temp):
        return self.convert_to_float(temp)


    def parse_mutation_event(self, mutation_event):
        return mutation_event[0], int(mutation_event[1:-1]), mutation_event[-1]


    def validate_mutation(self, pdb, mutation_event):
        # based on PoPMuSiC_2, Saraboji
        if pdb=="1LVE" and mutation_event=="V27BL":
            return mutation_event[0]+"29"+mutation_event[-1]
        elif pdb=="1LVE" and mutation_event=="L27CQ":
            return mutation_event[0]+"30"+mutation_event[-1]
        elif pdb == "1LVE" and mutation_event=="Y27DD":
            return mutation_event[0]+"31"+mutation_event[-1]
        # based on ThermoMuteDB_single
        if pdb=="1C9O" and mutation_event=="L67AL":
            return mutation_event[0]+"66"+mutation_event[-1]
        elif pdb=="1C9O" and mutation_event=="L66LA":
            return mutation_event[0]+"66"+mutation_event[-1]
        elif pdb=="1IMQ" and mutation_event=="T27TG":
            return mutation_event[0]+"27"+mutation_event[-1]
        elif pdb=="'-" and mutation_event=="AVIG":
            return mutation_event[0]+"6"+mutation_event[-1]
        elif pdb=="'-" and mutation_event=="AIIG":
            return mutation_event[0]+"2"+mutation_event[-1]
        else:
            return mutation_event


    def get_mutations(self, row):
        pass