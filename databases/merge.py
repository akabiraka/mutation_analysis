import sys
sys.path.append("../mutation_analysis")
import os
import pandas as pd

from objects.Mutation import Mutation

def append_non_redundent_mutation(inp_file_path):
    out_file_path = "data/merged.csv"

    df = pd.read_csv(inp_file_path)
    n_mutation_appened, n_mutation_updated = 0, 0
    for row in df.itertuples():
        mutation = Mutation().get_mutation(row)
        is_updated, is_appended = mutation.save(out_file_path)
        if is_appended: 
            n_mutation_appened += 1
            print("crnt_row_idx: {}, n_mutation_appened: {}".format(row.Index, n_mutation_appened))
            # break
        elif is_updated:
            n_mutation_updated += 1
            print("crnt_row_idx: {}, n_mutation_updated: {}".format(row.Index, n_mutation_updated))
            # break

    print("n_mutation_appened: {}".format(n_mutation_appened))
    print("n_mutation_updated: {}".format(n_mutation_updated))
    print("n_already_exists: {}".format(df.shape[0] - n_mutation_appened - n_mutation_updated))
append_non_redundent_mutation("data/clean_1/ProThermDB.csv")