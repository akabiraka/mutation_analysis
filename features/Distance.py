import sys
sys.path.append("../mutation_analysis")

from Bio.PDB import PDBParser
import numpy as np

class Distance(object):
    def __init__(self) -> None:
        pass

    def get(self, pdb_file:str, chain_id:str, central_residue_id:int, neighbor_residue_id:int, central_residue_atom="CB", neighbor_residue_atom="CB"):
        residues = PDBParser(QUIET=True).get_structure("", pdb_file)[0][chain_id]
        central_residue = residues[central_residue_id]
        neighbor_residue = residues[neighbor_residue_id]

        distance_vector = central_residue[central_residue_atom].coord - neighbor_residue[neighbor_residue_atom].coord
        distance = np.linalg.norm(distance_vector) #or np.sqrt(np.sum(distance_vector * distance_vector))
        unit_distance_vector = distance_vector / np.linalg.norm(distance_vector)

        return distance, distance_vector, unit_distance_vector



distance = Distance()
distance, distance_vector, unit_distance_vector = distance.get("data/pdbs_clean/1a0fA.pdb", "A", 5, 4)
print(distance, distance_vector, unit_distance_vector)
