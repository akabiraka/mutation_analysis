import sys
sys.path.append("../mutation_analysis")

import features.static as STATIC
from features.PSSM import PSSM
from encoding_methods.Onehot import Onehot
from encoding_methods.WaveEncoding import WaveEncoding

import numpy as np
from Bio import SeqIO
import itertools

class SAAFEC_SEQ(object):
    def __init__(self) -> None:
        self.pssm = PSSM()
        self.onehot = Onehot()
        self.waveenc = WaveEncoding()
        self.mutation_types = {"".join(x):i+1 for i, x in enumerate(itertools.permutations("ARNDCQEGHILKMFPSTWYV", 2))} #range(1, 380+1)
        # A:Aliphatic, B:aromatic, C:sulfur, D:hydroxyl, E:basic, F:acidic, G:amide
        self.chem_prop_dict = {"A":"A", "G":"A", "I":"A", "L":"A", "P":"A", "V":"A", 
                               "F":"B", "W":"B", "Y":"B", 
                               "C":"C", "M":"C",
                               "S":"D", "T":"D",
                               "R":"E", "H":"E", "K":"E", 
                               "D":"F", "E":"F", 
                               "N":"G", "Q":"G"}
        self.chem_prop_change_types = {"".join(x):i+1 for i, x in enumerate(itertools.product("ABCDEFG", "ABCDEFG"))} #range(1, 49+1)
        # A:small, B:medium, C:large
        self.size_dict = {"G":"A", "A":"A", "S":"A", "C":"A", "D":"A", "P":"A", "N":"A", "T":"A",
                          "Q":"B", "E":"B", "H":"B", "V":"B",   
                          "R":"C", "I":"C", "L":"C", "K":"C", "M":"C", "F":"C", "W":"C", "Y":"C"}
        self.size_change_types = {"".join(x):i+1 for i, x in enumerate(itertools.product("ABC", "ABC"))} #range(1, 9+1)
        # A:polar, B:nonpolar
        self.polarity_dict = {"R":"A", "N":"A", "D":"A", "Q":"A", "E":"A", "H":"A", "K":"A", "S":"A", "T":"A", "Y":"A", 
                              "A":"B", "C":"B", "G":"B", "I":"B", "L":"B", "M":"B", "F":"B", "P":"B", "W":"B", "V":"B"}
        self.polarity_change_types = {"".join(x):i+1 for i, x in enumerate(itertools.product("AB", "AB"))} #range(1, 4+1)
        # A:donor, B:acceptor, C:donor and acceptor, D:none
        self.hydrogen_bond_dict = {"R":"A", "K":"A", "W":"A", 
                                   "D":"B", "E":"B", 
                                   "N":"C", "Q":"C", "H":"C", "S":"C", "T":"C", "Y":"C",
                                   "A":"D", "C":"D", "G":"D", "I":"D", "L":"D", "M":"D", "F":"D", "P":"D", "V":"D"}
        self.hydrogen_bond_change_types = {"".join(x):i+1 for i, x in enumerate(itertools.product("ABCD", "ABCD"))} #range(1, 16+1)
        # A:hydrophobic, B:neutral, C:hydrophilic
        self.hydrophobicity_dict = {"A":"A", "C":"A", "I":"A", "L":"A", "M":"A", "F":"A", "W":"A", "V":"A", 
                                    "G":"B", "H":"B", "P":"B", "S":"B", "T":"B", "Y":"B", 
                                    "R":"C", "N":"C", "D":"C", "Q":"C", "E":"C", "K":"C"}
        self.hydrophobicity_change_types = {"".join(x):i+1 for i, x in enumerate(itertools.product("ABC", "ABC"))} #range(1, 9+1)

    def get_pseudo_PSSM(self, pssm_file, window_size=7):
        """See section 4.2.1 from paper.
        Args:
            pssm_file (_type_): _description_
            window_size (int, optional): 0 < window_size < L. Defaults to 7.
        Returns:
            ndarray: 1D. 20+20*window
        """
        df, _ = self.pssm.parse_pssm_output_file(pssm_file)
        P = np.array(df[range(2, 21+1)], dtype=np.float32) #shape: Lx20, P=log-odds
        # Equation 1
        P_bars = np.average(P, axis=0) #shape: 20, avg per column
        # Equation 2
        L = P.shape[0]
        pseudo_features = []
        for j in range(20):
            for k in range(1, window_size+1):
                x = 0.0
                for i in range(L-k): 
                    x += pow(P[i, j] - P[i+k, j], 2)
                x = x/(L-k)
                pseudo_features.append(x)
        features = np.hstack([P_bars, pseudo_features]) #shape: 20+20xwindow_size
        return features


    def get_neighbor_conservation_scores(self, pssm_file, i, n_neighbors=3):
        """See section 4.2.2 from paper.
        Args:
            pssm_file (str): file path
            i (int): zero_based_mutation_site
            n_neighbors (int, optional): Defaults to 7.
        Returns:
            ndarray: 1D. 20*n_neighbors
        """
        features = self.pssm.get_neighbor_logodds_by_mutation_site_in_middle(i=i, neighbors=n_neighbors, pssm_file=pssm_file) #shape: 20xwindow_size=140
        # print(features)
        features = np.hstack(features)
        return features
        

    def get_seq_neighbor_features(self, fasta_file, i, n_neighbors=3, smooth=True):
        """See section 4.2.3 from paper.
        Args:
            fasta_file (str): file path
            i (int): 0-based mutation site
            n_neighbors (int, optional): on both sides of i. Defaults to 3.
            smooth (bool, optional): The encoding type. Defaults to True.

        Returns:
            ndarray: [(2*n_neighbors+1)*20]
        """
        fasta = next(SeqIO.parse(fasta_file, "fasta"))
        seq = str(fasta.seq)
        features = self.onehot.mutation_site_with_neighbors(seq, i, n_neighbors, smooth)
        features = np.hstack(features)
        # print(features.shape)
        return features

    def property_change(self, wild_res, mut_res, prop_dict, change_types_dict):
        return change_types_dict[prop_dict[wild_res] + prop_dict[mut_res]]

    def get_physiochemical_properties_at_mutation_site(self, wild_residue, mutant_residue):
        net_vol = STATIC.RESIDUE_VOLUME[mutant_residue] - STATIC.RESIDUE_VOLUME[wild_residue]
        net_hydrophobicity = STATIC.MOON_HYDROPHOBIC_INDEX[mutant_residue] - STATIC.MOON_HYDROPHOBIC_INDEX[wild_residue]
        mutation_type = self.mutation_types[wild_residue + mutant_residue]
        # net_flexibility = ? # did not understand
        chem_prop_change = self.property_change(wild_residue, mutant_residue, self.chem_prop_dict, self.chem_prop_change_types)
        size_change = self.property_change(wild_residue, mutant_residue, self.size_dict, self.size_change_types)
        polarity_change = self.property_change(wild_residue, mutant_residue, self.polarity_dict, self.polarity_change_types)
        hydrogen_bond_change = self.property_change(wild_residue, mutant_residue, self.hydrogen_bond_dict, self.hydrogen_bond_change_types)
        hydrophobicity_change = self.property_change(wild_residue, mutant_residue, self.hydrophobicity_dict, self.hydrophobicity_change_types)
        features = np.array([net_vol, net_hydrophobicity, mutation_type, chem_prop_change, size_change, polarity_change, hydrogen_bond_change, hydrophobicity_change], dtype=np.float32)
        # print(features.shape)
        return features

    def extra(self, wild_residue, mutant_residue):
        continuous_mutation_type = self.waveenc.get(self.mutation_types[wild_residue+mutant_residue], dim=10)
        print(continuous_mutation_type)


saafec_seq = SAAFEC_SEQ()
# saafec_seq.get_pseudo_PSSM("data/pssms_s2648_wild/1a5eA.pssm")
# features = saafec_seq.get_neighbor_conservation_scores("data/pssms_s2648_wild/1a5eA.pssm", 0, 3)
# features = saafec_seq.get_seq_neighbor_features("data/fastas/1a5eA.fasta", 0, 3, True) # 1st amino acid is the mutation site
# features = saafec_seq.get_seq_neighbor_features("data/fastas/1a5eA.fasta", 155, 3, True) # 155th (last) amino acid is the mutation site
saafec_seq.get_physiochemical_properties_at_mutation_site("R", "A")
# saafec_seq.extra("R", "A")