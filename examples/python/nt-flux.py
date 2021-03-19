import sys
from math import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# import DiskModel class from disk_model.py
root_path = '../..'
sys.path.append(root_path)
from disk_model import DiskModel

bh_mass = 10.0
mdot_min = 0.01
mdot_max = 2.0


for spin in [0.0, 0.90, 0.99]:
    m = DiskModel(root_path+'/models/nt/disk-nt.so', bh_mass, spin, 'mdot=0.1')
    R = np.logspace(log10(m.r_min), log10(1e3), 100)
    F = np.zeros(R.shape)
    for i in range(len(R)): F[i] = m.flux(R[i])
    plt.loglog(R, F);
    m.close()
#end for


# plot flux
plt.xlabel(r'Radius [$r_{\rm g}$]')
plt.ylabel(r'Local flux [erg/cm$^2$/s]')
plt.show()

