import sys
sys.path.append("../mutation_analysis")

from objects.PDBData import PDBData
import pickle
  

def get_row_items(row):
    return row["pdb_id"].lower()[:4], row["chain_id"], row["mutation"], int(row["mutation_site"]), row["wild_residue"], row["mutant_residue"], row["ddG"]

def save_as_pickle(data, path):
    with open(path, "wb") as f: pickle.dump(data, f)

def load_pickle(path):
    with open(path, "rb") as f: return pickle.load(f)