import sys
sys.path.append("../mutation_analysis")

class ACDC_NN(object):
    def __init__(self) -> None:
        pass

    def get_variation(self, aa, is_wildtype=True):
        X = {"A": 0., "R": 0., "N": 0., "D": 0., "C": 0., "Q": 0., "E": 0., "G": 0., "H": 0., "I": 0.,
             "L": 0., "K": 0., "M": 0., "F": 0., "P": 0., "S": 0., "T": 0., "W": 0., "Y": 0., "V": 0.}
        if is_wildtype: X[aa] = -1.0
        else: X[aa] = 1.0
        return list(X.values())

acdc_nn = ACDC_NN()
v = acdc_nn.get_variation("W", True)
print(v)