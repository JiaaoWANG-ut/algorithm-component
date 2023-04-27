import numpy as np
from ase import Atoms, neighborlist
from ase.io import read

# 读取CIF文件
filename = "1.cif"  # 用你的CIF文件名替换
atoms = read(filename)

# 计算邻近原子列表
custom_cutoff_radius = 4.0  # 更改为所需的截断半径值
cutoff = np.full(len(atoms), custom_cutoff_radius)
#cutoff = neighborlist.natural_cutoffs(atoms)
nl = neighborlist.NeighborList(cutoff, self_interaction=False, bothways=True)
nl.update(atoms)

# 遍历所有原子，找到每个原子的近邻原子列表
for atom in atoms:
    atom_index = atom.index
    atom_symbol = atom.symbol
    atom_position = atom.position

    # 获取邻近原子索引
    neighbors_indices = nl.get_neighbors(atom_index)[0]

    print(f"原子 {atom_index} ({atom_symbol}):")
    print(f"  位置: {atom_position}")

    # 输出邻近原子的元素符号和位置
    print(f"  近邻原子:")
    for neighbor_index in neighbors_indices:
        neighbor_atom = atoms[neighbor_index]
        neighbor_symbol = neighbor_atom.symbol
        neighbor_position = neighbor_atom.position
        dist=np.linalg.norm(neighbor_position-atom_position)
        if dist < 3 and neighbor_symbol == 'O':
            print(f"    原子 {neighbor_index} ({neighbor_symbol}):")
            print(f"      位置: {neighbor_position}",dist)
