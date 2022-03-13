import sys
sys.path.append("../mutation_analysis")

from encoding_methods.WaveEncoding import WaveEncoding
import numpy as np
import matplotlib.pyplot as plt

Enc = WaveEncoding()
plt.figure(figsize=(20, 5))
n, dim = 100, 20
for i in range(n):
    y = np.array([Enc.get(i, dim) for i in range(100)])
    plt.plot(np.arange(n), y[:, 4:8])
    plt.scatter(np.repeat(i, dim), y[i])
    plt.legend(["dim %d"%p for p in range(4, 8)])
plt.show()