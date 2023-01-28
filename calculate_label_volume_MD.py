import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import numpy as np
import matplotlib.pyplot as plt

### settings: ###
lbl = 'ALX' # change to your label name
random_size = 10000 # benchmarked over a revalent range for a typical fluorescent label
cutoff = 3 # a label and a random point would match if the distance shorther than this cutoff
L = 30 # sets a box around the label's CA to generate random numbers

u = mda.Universe('./structure.pdb', './trajectory.xtc')

label_ca = u.select_atoms("resname " + lbl +" and name CA").positions[0]

label_xyz = []

for ts in u.trajectory[::10]:
    label = u.select_atoms("resname " + lbl)
    label_xyz.append(label.positions)

vols = []

random_x = np.random.uniform(label_ca[0] - L, label_ca[0] + L, random_size)
random_y = np.random.uniform(label_ca[1] - L, label_ca[1] + L, random_size)
random_z = np.random.uniform(label_ca[2] - L, label_ca[2] + L, random_size)
box_vol = (2*L)**3
radom_xyz = np.vstack([random_x, random_y, random_z]).T

for i in range(1, len(label_xyz)):
    dist_arr = np.zeros((np.vstack(label_xyz[:i]).shape[0], radom_xyz.shape[0]))
    distance_array(np.vstack(label_xyz[:i]), radom_xyz, result = dist_arr)
    count = ((dist_arr < cutoff).sum(axis = 0) > 0).sum()
    p = count/random_size
    vols.append(p * (2*L)**3)
    
plt.plot(vols, color = 'k', lw = 2)
plt.xlabel('Simulation time')
plt.ylabel(u'Volume, $Å^3$')
plt.show()
