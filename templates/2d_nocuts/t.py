import numpy as np
import json
lumi = 140
prefix = "../2d_cuts/temps_2D"
ts = []
dx = np.array([0.05, 2 * np.pi / 40])
print(dx)
for i in range(9):
    with open(prefix + str(i) + ".json", "r") as f:
        ts.append(json.load(f)["zvals"])       
ts = np.array(ts) * np.prod(dx) * lumi
print(np.sum(ts))
