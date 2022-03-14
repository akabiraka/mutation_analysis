import sys
sys.path.append("../mutation_analysis")
import pandas as pd


def print_features_statistics(inp_file, out_file):
        df = pd.read_csv(inp_file, header=None)
        stat_df = df.describe().transpose()
        stat_df.to_csv(out_file, index=False)


# print_features_statistics(inp_file="data/computed_features/SCONES_features_on_PoPMuSiC_2.csv", 
#                           out_file="data/computed_features/SCONES_features_on_PoPMuSiC_2_statistics.csv")

print_features_statistics(inp_file="data/computed_features/SAAFEC_SEQ_features_on_PoPMuSiC_2.csv", 
                          out_file="data/computed_features/SAAFEC_SEQ_features_on_PoPMuSiC_2_statistics.csv")                          
                                   