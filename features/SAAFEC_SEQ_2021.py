import sys
sys.path.append("../mutation_analysis")
from features.PSSM import PSSM
import numpy as np

class SAAFEC_SEQ(object):
    def __init__(self) -> None:
        self.pssm = PSSM()

    def compute_pseudo_PSSM(self, pssm_file, window_size=7):
        """_summary_

        Args:
            pssm_file (_type_): _description_
            window_size (int, optional): 0 < window_size < L. Defaults to 7.
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

saafec_seq = SAAFEC_SEQ()
saafec_seq.compute_pseudo_PSSM("data/pssms_s2648_wild/1a5eA.pssm")