import sys
sys.path.append("../mutation_analysis")

from Bio.PDB import PDBParser
import numpy as np

class Distance(object):
    def __init__(self) -> None:
        pass

    def get(self, pdb_file:str, chain_id:str, center_residue_id:int, neighbor_residue_id:int, central_residue_atom="CB", neighbor_residue_atom="CB"):
        chain = PDBParser(QUIET=True).get_structure("", pdb_file)[0][chain_id]
        central_residue = chain[center_residue_id]
        neighbor_residue = chain[neighbor_residue_id]

        if central_residue.resname=="GLY": central_residue_atom="CA"
        if neighbor_residue.resname=="GLY": neighbor_residue_atom="CA"

        if central_residue.has_id(central_residue_atom)==False:
            if central_residue_atom=="CB": central_residue_atom="CA"
            elif central_residue_atom=="CA": central_residue_atom="CB"
            else: raise NotImplementedError()

        if neighbor_residue.has_id(neighbor_residue_atom)==False:
            if neighbor_residue_atom=="CB": neighbor_residue_atom="CA"
            elif neighbor_residue_atom=="CA": neighbor_residue_atom="CB"
            else: raise NotImplementedError()

        distance_vector = central_residue[central_residue_atom].coord - neighbor_residue[neighbor_residue_atom].coord
        distance = np.linalg.norm(distance_vector) #or np.sqrt(np.sum(distance_vector * distance_vector))
        unit_distance_vector = distance_vector / np.linalg.norm(distance_vector)

        return distance, distance_vector, unit_distance_vector



# distance = Distance()
# distance, distance_vector, unit_distance_vector = distance.get("data/pdbs_clean/5dfrA.pdb", "A", (' ', 122, ' '), (' ', 123, ' '))
# print(distance, distance_vector, unit_distance_vector)
