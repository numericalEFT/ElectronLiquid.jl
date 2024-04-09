
import numpy as np
import matplotlib.pyplot as plt

mass2 = np.array([1e-2,1e-4,1e-6,1e-8])
val = np.array([1.179, 1.379, 1.422, 1.422])
err = np.array([0.003,0.005,0.007,0.008])

plt.figure()
plt.xscale("log")
plt.errorbar(mass2, val, yerr=err, fmt="o-", markersize=4)
plt.savefig("testPP.pdf")
