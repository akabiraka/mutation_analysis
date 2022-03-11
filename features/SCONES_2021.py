import sys
sys.path.append("../mutation_analysis")

from Bio.PDB.Polypeptide import three_to_one
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser

import features.static as STATIC
from features.SS_SA_RASA import SS_SA_RASA
from features.Distance import Distance
from features.Orientation import Orientation

class SCONES(object):
    def __init__(self, computed_features_file="data/computed_features/SCONES_features_on_PoPMuSiC_2.csv") -> None:
        self.ss_sa_rasa = SS_SA_RASA()
        self.distance = Distance()
        self.orentation = Orientation()
        self.computed_features_file = computed_features_file

    def print_features_statistics(self):
        df = pd.read_csv(self.computed_features_file, header=None)
        stat_df = df.describe().transpose()
        stat_df.to_csv("data/computed_features/SCONES_features_on_PoPMuSiC_2_statistics.csv", index=False)
        

    def get_all_features(self, cln_pdb_file, chain_id, center_residue_id, neighbor_residue_id, central_residue_atom="CB", neighbor_residue_atom="CB"):
        center_residue_features = self.get_node_features(cln_pdb_file, chain_id, center_residue_id)
        neighbor_residue_features = self.get_node_features(cln_pdb_file, chain_id, neighbor_residue_id)
        edge_features = self.edge_features(cln_pdb_file, chain_id, center_residue_id, neighbor_residue_id, central_residue_atom, neighbor_residue_atom)
        features = np.hstack([center_residue_features, neighbor_residue_features, edge_features])
        # print(features.shape, features)
        return features

    def get_node_features(self, cln_pdb_file, chain_id, residue_id):
        residue = PDBParser(QUIET=True).get_structure(" ", cln_pdb_file)[0][chain_id][residue_id]
        resname_1_letter = three_to_one(residue.resname)
        formal_charge = STATIC.AA_FORMAL_CHARGE[resname_1_letter]
        normalized_van_der_waals_vol = STATIC.NORMALIZED_VAN_DER_WAALS_VOL[resname_1_letter]
        hydropathy = STATIC.HYDROPATHY_INDEX[resname_1_letter]
        steric = STATIC.STERIC_PARAMETER[resname_1_letter]
        polarity = STATIC.POLARITY[resname_1_letter]
        rasa_tripeptide = STATIC.RASA_TRIPEPTIDE[resname_1_letter]
        _, sasa = self.ss_sa_rasa.get_ss_and_rasa_at_residue(cln_pdb_file, chain_id, residue_id)

        features = np.array([formal_charge, normalized_van_der_waals_vol, hydropathy, steric, polarity, rasa_tripeptide, sasa], dtype=np.float32)
        # print(features)
        return features

    def edge_features(self, cln_pdb_file, chain_id, center_residue_id, neighbor_residue_id, central_residue_atom="CB", neighbor_residue_atom="CB"):
        dist_features = self.get_distance_features(cln_pdb_file, chain_id, center_residue_id, neighbor_residue_id, central_residue_atom, neighbor_residue_atom)
        orientation_features = self.get_orientation_features(cln_pdb_file, chain_id, center_residue_id, neighbor_residue_id)
        features = np.hstack([dist_features, orientation_features])
        # print(features)
        return features

    def get_distance_features(self, cln_pdb_file, chain_id, center_residue_id, neighbor_residue_id, central_residue_atom="CB", neighbor_residue_atom="CB"):
        d, _, _ = self.distance.get(cln_pdb_file, chain_id, center_residue_id, neighbor_residue_id, central_residue_atom, neighbor_residue_atom)
        features = np.array([d, pow(d, -1), pow(d, -2), pow(d, -3), np.exp(-d/1), np.exp(-d/.8)], dtype=np.float32)
        # print(features)
        return features

    def get_orientation_features(self, cln_pdb_file, chain_id, center_residue_id, neighbor_residue_id):
        phi, psi, omega = self.orentation.get_dihedral_angles(cln_pdb_file, chain_id, center_residue_id, neighbor_residue_id)
        features = np.array([np.sin(phi), np.cos(phi), np.sin(psi), np.cos(psi), np.sin(omega), np.cos(omega)], dtype=np.float32)
        # print(features)
        return features
    
    


scones = SCONES()
scones.print_features_statistics()
# scones.get_distance_features("data/pdbs_clean/1lveA.pdb", "A", (" ", 5, " "), (" ", 4, " "), "CB", "CB")
# scones.get_orientation_features("data/pdbs_clean/1lveA.pdb", "A", (" ", 5, " "), (" ", 4, " "))
# scones.edge_features("data/pdbs_clean/1lveA.pdb", "A", (" ", 5, " "), (" ", 4, " "), "CB", "CB")
# scones.get_all_features("data/pdbs_clean/1lveA.pdb", "A", (" ", 5, " "), (" ", 4, " "), "CB", "CB")