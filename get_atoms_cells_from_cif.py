from ase.io import read
import numpy as np

atoms = read('2.cif')
cell = np.array(atoms.get_cell())
positions = atoms.get_positions()
symbols = atoms.get_chemical_symbols()
atoms = []
for symbol, position in zip(symbols, positions):
    atoms.append([symbol, position])

#cell
#atoms
