#!/usr/bin/python3

import sys
import numpy as np
import matplotlib.pyplot as plt

if "-h" in sys.argv:
    print("""Usage: ./prog [options]
        -h         shows help
        -b [ínt]   sets number of bars in the plot
    """)

breaks = 50
if "-b" in sys.argv:
    breaks = int(sys.argv[sys.argv.index("-b") + 1])

# Read data from pipe
data = sys.stdin.read().split()

# Create numpy array
data = np.array([float(i) for i in data if i is not "" ])

plt.hist(data, bins=breaks, color="purple")
plt.xlabel("Execution time (s)")
plt.ylabel("Absolute Frequency")
plt.savefig("hist.svg")
plt.show()
