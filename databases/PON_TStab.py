import sys
sys.path.append("../mutation_analysis")

import pandas as pd

class PON_TStab():
    def __init__(self) -> None:
        self.file_path = "data/downloaded_as/PON-TStab_dataset.xlsx"
        self.df = pd.read_excel(self.file_path)
        print(self.df.shape)


ponTStab = PON_TStab()