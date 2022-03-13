import sys
sys.path.append("../mutation_analysis")

import features.static as STATIC
from features.PSSM import PSSM
from encoding_methods.Onehot import Onehot

import numpy as np
from Bio import SeqIO
import itertools

class SAAFEC_SEQ(object):
    def __init__(self) -> None:
        self.pssm = PSSM()
        self.onehot = Onehot()
        self.mutation_types = {"".join(x):i+1 for i, x in enumerate(itertools.permutations("ARNDCQEGHILKMFPSTWYV", 2))} #range(1, 380+1)


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

    def get_physiochemical_properties_at_mutation_site(self, wild_residue, mutant_residue):
        net_vol = STATIC.RESIDUE_VOLUME[mutant_residue] - STATIC.RESIDUE_VOLUME[wild_residue]
        net_hydrophobicity = STATIC.MOON_HYDROPHOBIC_INDEX[mutant_residue] - STATIC.MOON_HYDROPHOBIC_INDEX[wild_residue]
        mutation_type = self.mutation_types[wild_residue+mutant_residue]
        print(net_vol, net_hydrophobicity, mutation_type)


saafec_seq = SAAFEC_SEQ()
# saafec_seq.get_pseudo_PSSM("data/pssms_s2648_wild/1a5eA.pssm")
# features = saafec_seq.get_neighbor_conservation_scores("data/pssms_s2648_wild/1a5eA.pssm", 0, 3)
# features = saafec_seq.get_seq_neighbor_features("data/fastas/1a5eA.fasta", 0, 3, True) # 1st amino acid is the mutation site
# features = saafec_seq.get_seq_neighbor_features("data/fastas/1a5eA.fasta", 155, 3, True) # 155th (last) amino acid is the mutation site
saafec_seq.get_physiochemical_properties_at_mutation_site("R", "A")