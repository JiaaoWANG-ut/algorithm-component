import sys
import numpy as np
from ase import Atoms, neighborlist
from ase.io import read
import re
import fnl 

import time
np.set_printoptions(threshold=sys.maxsize)

# Input file 
filename = "example.cif"  
atoms = read(filename)

# Cutoff radius 
rcut = 6.0  
symbols= re.split("\d+", atoms.get_chemical_formula())
nelement = len(list(filter(bool, symbols))) 
pos_car = atoms.get_positions()
cell = atoms.cell.array
atomicNrs = np.array(atoms.get_atomic_numbers())
uniqueNrs = []
for i in atomicNrs:
    if i not in uniqueNrs:
        uniqueNrs.append(i)
unique = np.array(uniqueNrs)

start = time.time()
fnl.nlist.update_nlist(atomicNrs, uniqueNrs, pos_car, cell, rcut)
end = time.time()
print ('Time for neighborlist calculation :', end - start)
print ('total neighbor list=', fnl.nlist.tneighs)
print ('total neighbor list in a cell=', fnl.nlist.tneighs_incell)
print ('total number of neighbors=', fnl.nlist.tnum_neigh)
print ('position of the total atoms (including ghost atoms)=', fnl.nlist.pool_pos_car)
