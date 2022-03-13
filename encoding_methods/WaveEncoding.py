import sys
sys.path.append("../mutation_analysis")
import numpy as np
import math

class WaveEncoding(object):
    def __init__(self) -> None:
        pass

    def get(self, i, dim):
        enc = np.zeros(dim)
        div_term = np.exp(np.arange(0, dim, 2) * -(math.log(10000.0) / dim))
        enc[0::2] = np.sin(np.array([i]) * div_term)
        enc[1::2] = np.cos(np.array([i]) * div_term)
        return enc
