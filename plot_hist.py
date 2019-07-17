#!/usr/bin/python3

import sys
import numpy as np
import matplotlib.pyplot as plt

def argvNext(arg):
    return sys.argv[sys.argv.index(arg) + 1]

if "-h" in sys.argv:
    print("""Usage: ./prog [options]
        -h           shows help
        -b [ínt]     sets number of bars in the plot
        -f [string]  sets the output filename
    """)

breaks = 50
if "-b" in sys.argv:
    breaks = int(argvNext("-b"))

fname = "hist.svg"
if "-f" in sys.argv:
    fname = argvNext("-f")

# Read data from pipe
data = sys.stdin.read().split()

# Create numpy array
data = np.array([float(i) for i in data if i is not "" ])

plt.hist(data, bins=breaks, color="purple")
plt.xlabel("Execution time (s)")
plt.ylabel("Absolute Frequency")
plt.savefig(fname)
print("Saved file as '{}'".format(fname))
plt.show()
