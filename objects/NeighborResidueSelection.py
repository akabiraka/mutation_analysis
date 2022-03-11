import sys
sys.path.append("../mutation_analysis")

import numpy as np
from Bio.PDB import PDBParser

class NeighborResidueSelection(object):
    def __init__(self) -> None:
        pass

    def get_all_residue_id_vs_distance(self, pdb_file, chain_id, center_residue_id, atom, sort=True):
        chain = PDBParser(QUIET=True).get_structure("", pdb_file)[0][chain_id]
        center_residue = chain[center_residue_id]

        if center_residue.has_id(atom)==False and atom=="CB": atom="CA"
        if center_residue.resname=="GLY": atom="CA"

        all_residue_id_vs_distance = []
        for i, residue in enumerate(chain.get_residues()):
            if residue.has_id("CA")==False or residue.has_id("N")==False or residue.has_id("C")==False: continue #for all amino acids
            if residue.resname=="GLY": atom="CA"
            if residue.has_id("CB")==False: atom="CA"
            if center_residue.id == residue.id: continue #not taking the residue itself as its neighbor
            
            diff_vector = center_residue[atom].coord - residue[atom].coord
            dist = np.linalg.norm(diff_vector)
            all_residue_id_vs_distance.append([residue.id, dist])

        if not sort: return all_residue_id_vs_distance

        sorted_all_residue_id_vs_distance = np.array(sorted(all_residue_id_vs_distance, key=lambda x: x[1]), dtype=object)
        return sorted_all_residue_id_vs_distance

    def get_k_nearest_neighbor_residue_ids(self, pdb_file, chain_id, center_residue_id, atom="CB", k=15):
        """
        Args:
            pdb_file (str): file path
            chain_id (str): "A"
            center_residue_id (tuple): i.e ('H_GLC', 100, 'A') or 100
            atom (str, optional): Defaults to "CB".
                Options: CA, CB, C, N
            k (int, optional): Defaults to 15.
        """
        all_residue_id_vs_distance = self.get_all_residue_id_vs_distance(pdb_file, chain_id, center_residue_id, atom)
        return all_residue_id_vs_distance[:k, 0].tolist()

    def get_distance_based_neighbor_residue_ids(self, pdb_file, chain_id, center_residue_id, atom="CB", distance=8.0):
        """
        Args:
            pdb_file (str): file path
            chain_id (str): string of len 1
            center_residue_id (tuple): i.e ('H_GLC', 100, 'A') or 100
            distance (float, optional): in Angstram. Defaults to 8.
        """
        all_residue_id_vs_distance = self.get_all_residue_id_vs_distance(pdb_file, chain_id, center_residue_id, atom)
        # print(all_residue_id_vs_distance)

        # getting neighbors whose distance are less than input distance
        neighbor_residue_ids = []
        for residue_id_vs_distance in all_residue_id_vs_distance:
            if residue_id_vs_distance[1]<=distance:
                neighbor_residue_ids.append(residue_id_vs_distance[0])
        # print(neighbor_residue_ids)
        return neighbor_residue_ids


# nrs = NeighborResidueSelection()
# neighbor_residue_ids = nrs.get_distance_based_neighbor_residue_ids("data/pdbs_clean/1fnaA.pdb", "A", (" ", 8, " "), "CB", 8.0)
# print(neighbor_residue_ids)
# neighbor_residue_ids = nrs.get_k_nearest_neighbor_residue_ids("data/pdbs_clean/1lveA.pdb", "A", (" ", 27, " "), "CB", 15)
# print(neighbor_residue_ids)