import sys
sys.path.append("../mutation_analysis")
from Bio.PDB import PDBParser
from Bio.PDB.vectors import calc_dihedral

class Orientation(object):
    def __init__(self) -> None:
        pass

    def __get_residues(self, pdb_file, chain_id, center_residue_id, neighbor_residue_id):
        residues = PDBParser(QUIET=True).get_structure("", pdb_file)[0][chain_id]
        central_residue = residues[center_residue_id]
        neighbor_residue = residues[neighbor_residue_id]
        return central_residue, neighbor_residue


    def get_dihedral_angles(self, pdb_file, chain_id, center_residue_id, neighbor_residue_id):
        """used the following two links to compute phi, psi and omega angels
        https://biopython.org/docs/dev/api/Bio.PDB.internal_coords.html#Bio.PDB.internal_coords.IC_Residue.pick_angle
        https://biopython.org/docs/latest/api/Bio.PDB.vectors.html?highlight=calc_dihedral#Bio.PDB.vectors.calc_dihedral
        (sC, nN, nCA, nC)   # phi i+1
        (sN, sCA, sC, nN)   # psi
        (sCA, sC, nN, nCA)  # omega i+1
        The angles are in [-pi, pi]
        """

        central_residue, neighbor_residue = self.__get_residues(pdb_file, chain_id, center_residue_id, neighbor_residue_id)

        sN, sCA, sC = central_residue["N"].get_vector(), central_residue["CA"].get_vector(), central_residue["C"].get_vector()
        nN, nCA, nC = neighbor_residue["N"].get_vector(), neighbor_residue["CA"].get_vector(), neighbor_residue["C"].get_vector()
        
        phi = calc_dihedral(sC, nN, nCA, nC)
        psi = calc_dihedral(sN, sCA, sC, nN)
        omega = calc_dihedral(sCA, sC, nN, nCA)

        return phi, psi, omega


    def get_planar_angle(self, pdb_file, chain_id, central_residue_id, neighbor_residue_id):
        """Uses SCONES, 2021 Figure 4.
        """
        central_residue, neighbor_residue = self.__get_residues(pdb_file, chain_id, central_residue_id, neighbor_residue_id)
        
        sCA, sCB = central_residue["CA"].get_vector(), central_residue["CB"].get_vector()
        nCA, nCB = neighbor_residue["CA"].get_vector(), neighbor_residue["CB"].get_vector()
        
        planer = calc_dihedral(sCA, sCB, nCB, nCA)
        return planer


# orientation = Orientation()
# phi, psi, omega = orientation.get_dihedral_angles("data/pdbs_clean/1a0fA.pdb", "A", 5, 10)
# print(phi, psi, omega)

# phi, psi, omega = orientation.get_dihedral_angles("data/pdbs_clean/1a0fA.pdb", "A", 10, 5)
# print(phi, psi, omega)

# planer = orientation.get_planar_angle("data/pdbs_clean/1a0fA.pdb", "A", 10, 5)
# print(planer)
# planer = orientation.get_planar_angle("data/pdbs_clean/1a0fA.pdb", "A", 5, 10)
# print(planer)